// Created by Petr Karnakov on 25.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <cassert>
#include <map>
#include <unordered_map>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "distr.h"
#include "distr_particles.h"
#include "dump/dumper.h"
#include "util/format.h"
#include "util/mpi.h"

#include "CubismNoCopy/BlockInfo.h"
#include "CubismNoCopy/BlockLab.h"
#include "CubismNoCopy/BlockLabMPI.h"
#include "CubismNoCopy/Grid.h"
#include "CubismNoCopy/GridMPI.h"

#if USEFLAG(HDF)
#include "CubismNoCopy/HDF5Dumper_MPI.h"
#endif

#include "CubismNoCopy/StencilInfo.h"

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

// Hide implementation and avoid collision with GBlk
// TODO: rename cubism_impl::GBlk
namespace cubismnc_impl {

// Static parameters for Cubism
// bx, by, bz - block size
// halo - number of ghost cells uniform in all directions
template <class Scal_, size_t bx_, size_t by_, size_t bz_, int halo_>
struct GPar {
  using Scal = Scal_;
  static const size_t bx = bx_;
  static const size_t by = by_;
  static const size_t bz = bz_;
  static const int halo = halo_;
};

// Field wrapper to emulate Cubism block behavior
// Par_ - instance of GPar
template <class Par_, size_t Pad_ = 1>
struct GFieldViewRaw {
  using Par = Par_;
  using Scal = typename Par::Scal;

  // inner block strides (w/o halos)
  static const size_t bx = Par::bx;
  static const size_t by = Par::by;
  static const size_t bz = Par::bz;

  // required by Cubism
  static const int sizeX = bx;
  static const int sizeY = by;
  static const int sizeZ = bz;

  // strides including halos
  // XXX: [fabianw@mavt.ethz.ch; 2019-11-13] The +Pad_ is needed for efficient
  // transformation between indexing types (cell, face, node).
  static const size_t stridex = Par::bx + 2 * Par::halo + Pad_;
  static const size_t stridey = Par::by + 2 * Par::halo + Pad_;
  static const size_t stridez = Par::bz + 2 * Par::halo + Pad_;
  static const int halo =
      Par::halo; // assumes # of halos is the same in all dimensions

  // carried element type must be Scal, do not change this
  using Elem = Scal;
  using ElementType = Elem;
  using element_type = Elem;

  GFieldViewRaw() = delete;
  GFieldViewRaw(Elem* const base, const int idx[3], const int n_components = 1)
      : data(
            base +
            n_components * (halo + halo * stridex + halo * stridex * stridey))
      , data_halo(base)
      , index{idx[0], idx[1], idx[2]}
      , n_comp(n_components) {
    assert(data_halo != data); // halo = 0 will not work in this model
  }

  // SoA field
  Elem* const data; // block data
  Elem* const data_halo; // block data including halos

  const int index[3]; // corresponding block index

  // Number of components carried by field.  Parameter is not passed as
  // template parameter such that scalar and vector (field) views can be
  // combined in the same field compound (to maximize MPI communication buffer
  // size).
  const int n_comp;

  // block access emulator
  inline Elem& operator()(
      unsigned int x, unsigned int y = 0, unsigned int z = 0) {
    assert(data != nullptr);
    assert(x < (int)bx);
    assert(y < (int)by);
    assert(z < (int)bz);
    return data[n_comp * (x + stridex * (y + stridey * z))];
  }

  inline Elem& LinAccess(unsigned int x) {
    assert(data_halo != nullptr);
    assert(x < n_comp * stridex * stridey * stridez);

    return data_halo[n_comp * x];
  }
};

// TView - instance of GFieldViewRaw
template <class TView>
class LabPer : public BlockLab<TView> {
  using FieldView = TView;
  using ElementTypeBlock = typename FieldView::Elem;

 public:
  virtual inline std::string name() const {
    return "LabPer";
  }
  bool is_xperiodic() {
    return true;
  }
  bool is_yperiodic() {
    return true;
  }
  bool is_zperiodic() {
    return true;
  }
};

// TView - instance of GFieldViewRaw
template <class TView>
using GLab = BlockLabMPI<LabPer<TView>>;

// TView - instance of GFieldViewRaw
template <class TView>
using GGrid = GridMPI<Grid<TView>>;

// Par - instance of GPar
template <class Par, class M_>
class Cubismnc : public DistrMesh<M_> {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using RedOp = typename M::Op;

  Cubismnc(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var);
  ~Cubismnc();

 private:
  using FieldView = GFieldViewRaw<Par>;
  using Lab = GLab<FieldView>;
  using Grid = GGrid<FieldView>;
  using Elem = typename FieldView::Elem;
  using Synch = typename Grid::Synch;

  using P = DistrMesh<M>; // parent
  using MIdx = typename M::MIdx;
  using BlockInfoProxy = generic::BlockInfoProxy<M::dim>;

  // Converts Cubism BlockInfo to BlockInfoProxy
  static std::vector<BlockInfoProxy> Convert(
      const std::vector<BlockInfo>&, MIdx bs, size_t hl);

  std::vector<size_t> TransferHalos() override;
  std::vector<size_t> TransferHalos(bool inner) override;
  void ReduceSingleRequest(const std::vector<RedOp*>& blocks) override;
  void Bcast(const std::vector<size_t>& bb) override;
  void Scatter(const std::vector<size_t>& bb) override;
  void DumpWrite(const std::vector<size_t>& bb) override;
  void TransferParticles(const std::vector<size_t>& bb) override;

  using P::comm_;
  using P::dim;
  using P::domain_;
  using P::frame_;
  using P::isroot_;
  using P::kernelfactory_;
  using P::kernels_;
  using P::mshared_;
  using P::rank_from_id_;
  using P::stage_;
  using P::var;

  struct CheckProcs {
    CheckProcs(MPI_Comm comm, MIdx nprocs) {
      int commsize;
      MPI_Comm_size(comm, &commsize);
      fassert(
          commsize == nprocs.prod(), //
          util::Format(
              "Number of MPI tasks {} does not match the number of subdomains "
              "{}",
              commsize, nprocs));
    }
  };

  CheckProcs checkprocs_; // Enforces the check before construction of grid_
  Grid grid_;
  std::vector<std::vector<FieldView>> fviews_; // fields from last communication
  Synch* sync_;
  size_t n_fields_;
  std::map<MIdx, size_t> midx_to_kernel_;
  int commsize_;
};

// B_ - instance of GFieldViewRaw
template <class B_>
struct StreamHdfScal {
  using B = B_;
  using Scal = typename B::Scal;

  // Required by Cubism
  static std::string NAME;
  static std::string EXT;
  static const int NCHANNELS = 1;
  static const int CLASS = 0;

  B& b;
  const int aos;

  StreamHdfScal(B& b, const int idx) : b(b), aos(idx) {
    assert(aos >= 0);
  }

  // write
  void operate(const int ix, const int iy, const int iz, Scal* out) const {
    out[0] = (&b(ix, iy, iz))[aos];
  }

  // read
  void operate(const Scal* out, const int ix, const int iy, const int iz) {
    (&b(ix, iy, iz))[aos] = out[0];
  }

  static const char* getAttributeName() {
    return "Scalar";
  }
};

// B_ - instance of GFieldViewRaw
template <class B_>
struct StreamHdfVect {
  using B = B_;
  using Scal = typename B::Scal;

  // Required by Cubism
  static std::string NAME;
  static std::string EXT;
  static const int NCHANNELS = 3;
  static const int CLASS = 0;

  B& b;
  const int aos;

  StreamHdfVect(B& b, const int idx) : b(b), aos(idx) {
    assert(3 == b.n_comp);
    assert(0 == aos);
  }

  // write
  void operate(const int ix, const int iy, const int iz, Scal* out) const {
    const Scal* v = &((const Scal*)&b(ix, iy, iz))[aos];
    out[0] = v[0];
    out[1] = v[1];
    out[2] = v[2];
  }

  // read
  void operate(const Scal* out, const int ix, const int iy, const int iz) {
    Scal* v = &((Scal*)&b(ix, iy, iz))[aos];
    v[0] = out[0];
    v[1] = out[1];
    v[2] = out[2];
  }

  static const char* getAttributeName() {
    return "Vector";
  }
};

template <class Par, class M>
auto Cubismnc<Par, M>::Convert(
    const std::vector<BlockInfo>& infos, MIdx blocksize, size_t halos)
    -> std::vector<BlockInfoProxy> {
  std::vector<BlockInfoProxy> proxies;
  for (size_t i = 0; i < infos.size(); i++) {
    const BlockInfo& info = infos[i];
    BlockInfoProxy p;
    p.index = MIdx(info.index);
    p.cellsize = Vect(info.h_gridpoint);
    p.blocksize = blocksize;
    p.halos = halos;
    proxies.push_back(p);
  }
  return proxies;
}

template <class Par, class M>
Cubismnc<Par, M>::Cubismnc(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var)
    : DistrMesh<M>(comm, kf, var)
    , checkprocs_(comm, domain_.nprocs)
    , grid_(
          domain_.nprocs[0], domain_.nprocs[1], domain_.nprocs[2], //
          domain_.nblocks[0], domain_.nblocks[1], domain_.nblocks[2], //
          domain_.extent, comm) {
  fassert_equal(domain_.blocksize[0], FieldView::bx);
  fassert_equal(domain_.blocksize[1], FieldView::by);
  fassert_equal(domain_.blocksize[2], FieldView::bz);

  {
    MPI_Comm_size(comm, &commsize_);
    int rank;
    MPI_Comm_rank(comm, &rank);
    isroot_ = (0 == rank);
  }

  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Get rid of BlockInfo type
  const std::vector<BlockInfo> infos = grid_.getBlocksInfo();

  std::vector<BlockInfoProxy> proxies =
      Convert(infos, domain_.blocksize, domain_.halos);
  for (size_t i = 0; i < proxies.size(); ++i) {
    auto& p = proxies[i];
    p.isroot = (p.index == MIdx(0));
    p.islead = (i == 0);
    p.globalsize = domain_.nprocs * domain_.nblocks * domain_.blocksize;
    midx_to_kernel_[p.index] = i;
  }

  comm_ = grid_.getCartComm();
  this->MakeKernels(proxies);

  rank_from_id_ = [domain = domain_](int id) -> int {
    const MIdx global_blocks = domain.nblocks * domain.nprocs;
    const MIdx block = GIndex<int, dim>(global_blocks).GetMIdx(id);
    const MIdx proc = block / domain.nblocks;
    const MIdx proc_rev(proc[2], proc[1], proc[0]);
    const MIdx nprocs = domain.nprocs;
    const MIdx nprocs_rev(nprocs[2], nprocs[1], nprocs[0]);
    return GIndex<int, dim>(nprocs_rev).GetIdx(proc_rev);
  };
  for (auto& kernel : kernels_) {
    auto& m = kernel->GetMesh();
    m.SetHandlerMpiRankFromId(rank_from_id_);
  }
}

template <class Par, class M>
Cubismnc<Par, M>::~Cubismnc() {}

template <class Par, class M>
auto Cubismnc<Par, M>::TransferHalos() -> std::vector<size_t> {
  const std::vector<BlockInfo> infos = grid_.getBlocksInfo(); // all blocks
  fassert_equal(infos.size(), kernels_.size());
  auto& m = kernels_.front()->GetMesh();
  // Perform communication if necessary
  if (m.GetComm().size() > 0) {
    // 0. Construct vector of fields to communicate for all blocks on this rank
    std::vector<std::vector<FieldView>> fviews;
    const size_t n_fields = m.GetComm().size();
    fviews.resize(n_fields);
    for (size_t i = 0; i < infos.size(); ++i) {
      auto& reqs = kernels_[i]->GetMesh().GetComm();
      assert(n_fields == fviews.size());
      for (size_t fi = 0; fi < n_fields; ++fi) {
        auto& req = reqs[fi];
        fviews[fi].push_back(
            FieldView(req->GetBasePtr(), infos[i].index, req->GetSize()));
      }
    }

    // 1. Exchange halos in buffer mesh.
    // stencil type
    const bool is_tensorial = true;
    const int nhalo_start[3] = {domain_.halos, domain_.halos, domain_.halos};
    const int nhalo_end[3] = {domain_.halos, domain_.halos, domain_.halos};

    // schedule asynchronous communication
    Synch& s = grid_.sync(
        fviews, nhalo_start, nhalo_end, is_tensorial,
        var.Int["mpi_compress_msg"], var.Int["histogram"]);

    // Get all blocks synchronized
    const std::vector<BlockInfo> avail = s.avail_halo();
    fassert_equal(avail.size(), infos.size());

    // 2. Load exchanged halos into the local fields
    // FIXME: [fabianw@mavt.ethz.ch; 2019-12-07] This parallel structure is
    // likely less efficient for single threaded execution
#pragma omp parallel
    {
      Lab l;
      l.prepare(grid_, s);
#pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < infos.size(); ++i) {
        for (auto& field : fviews) {
          FieldView& block = field[i];
          l.load(block, field);
        }
      }
    }
  }

  std::vector<size_t> bb(infos.size());
  std::iota(bb.begin(), bb.end(), 0);
  return bb;
}

template <class Par, class M>
auto Cubismnc<Par, M>::TransferHalos(bool inner) -> std::vector<size_t> {
  const std::vector<BlockInfo> infos = grid_.getBlocksInfo(); // all blocks
  fassert_equal(infos.size(), kernels_.size());

  fassert_equal(
      mshared_->GetComm().size(), 0,
      ". Comm() on shared mesh is not implemented");

  std::vector<size_t> bb;
  if (inner) {
    n_fields_ = kernels_.front()->GetMesh().GetComm().size();
    // Perform communication if necessary
    if (n_fields_ > 0) {
      // 0. Construct vector of fields to communicate for all blocks on this
      // rank
      fviews_.clear();
      fviews_.resize(n_fields_);
      for (size_t i = 0; i < infos.size(); ++i) {
        auto& reqs = kernels_[i]->GetMesh().GetComm();
        assert(n_fields_ == reqs.size());
        for (size_t fi = 0; fi < n_fields_; ++fi) {
          auto& req = reqs[fi];
          fviews_[fi].push_back(
              FieldView(req->GetBasePtr(), infos[i].index, req->GetSize()));
        }
      }

      // 1. Exchange halos in buffer mesh.
      // stencil type
      const bool is_tensorial = true;
      const int nhalo_start[3] = {domain_.halos, domain_.halos, domain_.halos};
      const int nhalo_end[3] = {domain_.halos, domain_.halos, domain_.halos};

      // schedule asynchronous communication
      sync_ = &grid_.sync(
          fviews_, nhalo_start, nhalo_end, is_tensorial,
          var.Int["mpi_compress_msg"], var.Int["histogram"]);
      const std::vector<BlockInfo> avail = sync_->avail_inner();

      // Create vector of indices and save block info to map
      std::set<size_t> availset;
      for (const auto& info : avail) {
        const size_t i = midx_to_kernel_.at(MIdx(info.index));
        bb.push_back(i);
        availset.insert(i);
      }

      // 2. Load exchanged halos into the local fields
      // FIXME: [fabianw@mavt.ethz.ch; 2019-12-07] This parallel structure is
      // likely less efficient for single threaded execution. Imbalanced loop
      // due to inner-most if statement.  Unless the number of blocks per field
      // is not small, parallelizing over the field container is not a good
      // strategy as it contains only few fields
#pragma omp parallel
      {
        Lab l;
        l.prepare(grid_, *sync_);
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < infos.size(); ++i) {
          for (auto& field : fviews_) {
            FieldView& block = field[i];
            if (availset.count(i)) {
              l.load(block, field);
            }
          }
        }
      }
    } else { // no communication, return all blocks
      bb.resize(infos.size());
      std::iota(bb.begin(), bb.end(), 0);
    }
  } else { // not inner
    // Perform communication if necessary
    if (n_fields_ > 0) {
      const std::vector<BlockInfo> avail = sync_->avail_halo();
      std::set<size_t> availset;
      for (auto info : avail) {
        const size_t i = midx_to_kernel_.at(MIdx(info.index));
        bb.push_back(i);
        availset.insert(i);
      }

      // Load exchanged halos into the local fields
      // FIXME: [fabianw@mavt.ethz.ch; 2019-12-07] This parallel structure is
      // likely less efficient for single threaded execution. Imbalanced loop
      // due to inner-most if statement.  Unless the number of blocks per field
      // is not small, parallelizing over the field container is not a good
      // strategy as it contains only few fields
#pragma omp parallel
      {
        Lab l;
        l.prepare(grid_, *sync_);
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < infos.size(); ++i) {
          for (auto& field : fviews_) {
            FieldView& block = field[i];
            if (availset.count(i)) {
              l.load(block, field);
            }
          }
        }
      }
    } else { // no communication, return no blocks
    }
  }

  return bb;
}

template <class Par, class M>
void Cubismnc<Par, M>::Bcast(const std::vector<size_t>& bb) {
  using OpConcat = typename UReduce<Scal>::OpCat;
  auto& vfirst = kernels_.front()->GetMesh().GetBcast();

  if (!vfirst.size()) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetBcast().size(), vfirst.size());
  }

  for (size_t i = 0; i < vfirst.size(); ++i) {
    if (OpConcat* o = dynamic_cast<OpConcat*>(vfirst[i].get())) {
      std::vector<char> buf = o->Neutral();

      if (isroot_) {
        // Read from root block
        for (auto b : bb) {
          const auto& m = kernels_[b]->GetMesh();
          if (m.IsRoot()) {
            const OpConcat* ob = static_cast<OpConcat*>(m.GetBcast()[i].get());
            ob->Append(buf);
          }
        }
      }

      int size = buf.size();
      MPI_Bcast(&size, 1, MPI_INT, 0, comm_);
      buf.resize(size);
      MPI_Bcast(buf.data(), buf.size(), MPI_CHAR, 0, comm_);

      // Write to all blocks
      for (auto b : bb) {
        const auto& m = kernels_[b]->GetMesh();
        OpConcat* ob = static_cast<OpConcat*>(m.GetBcast()[i].get());
        ob->Set(buf);
      }
    } else {
      fassert(false, "Bcast: Unknown M::Op instance");
    }
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearBcast();
  }
}
template <class Par, class M>
void Cubismnc<Par, M>::Scatter(const std::vector<size_t>& bb) {
  auto& vfirst = kernels_.front()->GetMesh().GetScatter();

  if (!vfirst.size()) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetScatter().size(), vfirst.size());
  }

  for (size_t q = 0; q < vfirst.size(); ++q) {
    int recvcount;
    int sizes_recvcount;
    std::vector<Scal> rbuf;
    std::vector<int> sizes_rbuf;

    const MPI_Datatype mscal = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    if (isroot_) {
      std::vector<Scal> buf;
      std::vector<int> dis(commsize_, 0);
      std::vector<int> cnt(commsize_, 0);
      std::vector<int> sizes_buf;
      std::vector<int> sizes_dis(commsize_, 0);
      std::vector<int> sizes_cnt(commsize_, 0);
      // find root block
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        if (m.IsRoot()) {
          auto& req = m.GetScatter()[q];
          size_t i = 0;
          // concatenate data for all blocks in buf
          for (int rank = 0; rank < commsize_; ++rank) {
            dis[rank] = buf.size();
            sizes_dis[rank] = sizes_buf.size();
            // XXX assuming the same number of blocks on all ranks
            for (size_t k = 0; k < bb.size(); ++k) {
              auto& v = (*req.first)[i];
              buf.insert(buf.end(), v.begin(), v.end());
              sizes_buf.push_back(v.size());
              ++i;
            }
            sizes_cnt[rank] = sizes_buf.size() - sizes_dis[rank];
            cnt[rank] = buf.size() - dis[rank];
          }
        }
      }
      // data recvcount
      MPI_Scatter(cnt.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, 0, comm_);
      rbuf.resize(recvcount);
      // data
      MPI_Scatterv(
          buf.data(), cnt.data(), dis.data(), mscal, rbuf.data(), recvcount,
          mscal, 0, comm_);
      // sizes recvcount
      MPI_Scatter(
          sizes_cnt.data(), 1, MPI_INT, &sizes_recvcount, 1, MPI_INT, 0, comm_);
      sizes_rbuf.resize(sizes_recvcount);
      // sizes
      MPI_Scatterv(
          sizes_buf.data(), sizes_cnt.data(), sizes_dis.data(), MPI_INT,
          sizes_rbuf.data(), sizes_recvcount, MPI_INT, 0, comm_);
    } else {
      // data recvcount
      MPI_Scatter(nullptr, 0, MPI_INT, &recvcount, 1, MPI_INT, 0, comm_);
      rbuf.resize(recvcount);
      // data
      MPI_Scatterv(
          nullptr, nullptr, nullptr, mscal, rbuf.data(), recvcount, mscal, 0,
          comm_);
      // sizes recvcount
      MPI_Scatter(nullptr, 0, MPI_INT, &sizes_recvcount, 1, MPI_INT, 0, comm_);
      sizes_rbuf.resize(sizes_recvcount);
      // sizes
      MPI_Scatterv(
          nullptr, nullptr, nullptr, MPI_INT, sizes_rbuf.data(),
          sizes_recvcount, MPI_INT, 0, comm_);
    }

    // write to blocks on current rank
    size_t off = 0;
    for (size_t k = 0; k < bb.size(); ++k) {
      auto& v = *kernels_[bb[k]]->GetMesh().GetScatter()[q].second;
      v = std::vector<Scal>(
          rbuf.data() + off, rbuf.data() + off + sizes_rbuf[k]);
      off += sizes_rbuf[k];
    }
  }

  // Clear requests
  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearScatter();
  }
}

template <class Par, class M>
void Cubismnc<Par, M>::ReduceSingleRequest(const std::vector<RedOp*>& blocks) {
  using OpScal = typename UReduce<Scal>::OpS;
  using OpScalInt = typename UReduce<Scal>::OpSI;
  using OpConcat = typename UReduce<Scal>::OpCat;

  auto* firstbase = blocks.front();

  if (auto* first = dynamic_cast<OpScal*>(firstbase)) {
    auto buf = first->Neutral(); // result

    // Reduce over blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase);
      other->Append(buf);
    }

    MPI_Op mpiop;
    if (dynamic_cast<typename UReduce<Scal>::OpSum*>(first)) {
      mpiop = MPI_SUM;
    } else if (dynamic_cast<typename UReduce<Scal>::OpProd*>(first)) {
      mpiop = MPI_PROD;
    } else if (dynamic_cast<typename UReduce<Scal>::OpMax*>(first)) {
      mpiop = MPI_MAX;
    } else if (dynamic_cast<typename UReduce<Scal>::OpMin*>(first)) {
      mpiop = MPI_MIN;
    } else {
      fassert(false, "Unknown reduction");
    }
    MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    // Reduce over ranks
    MPI_Allreduce(MPI_IN_PLACE, &buf, 1, mt, mpiop, comm_);

    // Write results to all blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase);
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpScalInt*>(firstbase)) {
    auto buf = first->Neutral(); // result

    // Reduce over blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase);
      other->Append(buf);
    }

    MPI_Op mpiop;
    if (dynamic_cast<typename UReduce<Scal>::OpMinloc*>(first)) {
      mpiop = MPI_MINLOC;
    } else if (dynamic_cast<typename UReduce<Scal>::OpMaxloc*>(first)) {
      mpiop = MPI_MAXLOC;
    } else {
      fassert(false, "Unknown reduction");
    }

    MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE_INT : MPI_FLOAT_INT);

    // Reduce over all ranks
    MPI_Allreduce(MPI_IN_PLACE, &buf, 1, mt, mpiop, comm_);

    // Write results to all blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase);
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpConcat*>(firstbase)) {
    auto buf = first->Neutral();

    // Reduce over local blocks
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpConcat*>(otherbase);
      other->Append(buf);
    }

    int bufsize = buf.size();

    if (isroot_) {
      std::vector<int> sizes(commsize_);

      MPI_Gather(&bufsize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, comm_);

      int size_all = 0;
      std::vector<int> offsets = {0};
      for (auto& q : sizes) {
        size_all += q;
        offsets.push_back(offsets.back() + q);
      }
      offsets.pop_back();
      assert(sizes.size() == offsets.size());

      std::vector<char> buf_all(size_all);

      MPI_Gatherv(
          buf.data(), buf.size(), MPI_CHAR, buf_all.data(), sizes.data(),
          offsets.data(), MPI_CHAR, 0, comm_);

      // Write results to root block only (assume first)
      // FIXME get IsRoot() from kernel_ after using std::vector for kernels
      first->Set(buf_all);
    } else {
      MPI_Gather(&bufsize, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm_);

      MPI_Gatherv(
          buf.data(), buf.size(), MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR,
          0, comm_);
    }
    return;
  }
  fassert(false, "Unknown M::Op implementation");
}

template <class Par, class M>
void Cubismnc<Par, M>::DumpWrite(const std::vector<size_t>& bb) {
  auto& mfirst = kernels_.front()->GetMesh();
  if (mfirst.GetDump().size()) {
    std::string dumpformat = var.String["dumpformat"];
#if USEFLAG(HDF)
    if (dumpformat == "default") {
      dumpformat = "hdf";
    }
#endif

    if (dumpformat == "hdf") {
#if USEFLAG(HDF)
      // Create FieldView's for dump
      const size_t n_fields = mfirst.GetDump().size();
      std::vector<BlockInfo> infos = grid_.getBlocksInfo(); // all blocks
      std::vector<std::vector<FieldView>> fviews(n_fields);
      for (size_t i = 0; i < infos.size(); ++i) {
        auto& reqs = kernels_[i]->GetMesh().GetDump();
        assert(n_fields == reqs.size());
        for (size_t fi = 0; fi < n_fields; ++fi) {
          auto& req = reqs[fi].first;
          fviews[fi].push_back(
              FieldView(req->GetBasePtr(), infos[i].index, req->GetStride()));
        }
      }

      // Write dump
      auto& reqs = mfirst.GetDump();
      fassert_equal(reqs.size(), fviews.size());
      for (size_t fi = 0; fi < fviews.size(); ++fi) {
        const auto name = GetDumpName(reqs[fi].second, "", frame_);
        const int aos_idx = (reqs[fi].first)->GetIndex();
        auto& blocks = fviews[fi];
        if (0 <= aos_idx) {
          StreamHdfScal<FieldView>::NAME = reqs[fi].second;
          DumpHDF5_MPI<StreamHdfScal<FieldView>>(
              blocks, aos_idx, grid_, frame_, frame_, name, ".", Vect(0),
              mfirst.GetCellSize(), true);
        } else if (3 == blocks[0].n_comp) {
          StreamHdfVect<FieldView>::NAME = reqs[fi].second;
          DumpHDF5_MPI<StreamHdfVect<FieldView>>(
              blocks, aos_idx, grid_, frame_, frame_, name, ".", Vect(0),
              mfirst.GetCellSize(), true);
        } else {
          fassert(false, "DumpWrite(): Support only size 1 and 3");
        }
      }
      if (isroot_) {
        std::cerr << "Dump " << frame_ << ": format=" << dumpformat
                  << std::endl;
      }
      ++frame_;
#else
      fassert(false, "Compiled without HDF5. Use 'set string dumpformat raw'.");
#endif
    } else {
      P::DumpWrite(bb);
    }
  }
}

template <class Par, class M>
void Cubismnc<Par, M>::TransferParticles(const std::vector<size_t>&) {
  ::TransferParticles<M>(kernels_, comm_, rank_from_id_);
}

template <class B>
std::string StreamHdfScal<B>::NAME = "alpha";
template <class B>
std::string StreamHdfScal<B>::EXT = "";

template <class B>
std::string StreamHdfVect<B>::NAME = "alpha";
template <class B>
std::string StreamHdfVect<B>::EXT = "";

} // namespace cubismnc_impl

template <class Scal, size_t bx, size_t by, size_t bz, size_t hl>
using GPar = cubismnc_impl::GPar<Scal, bx, by, bz, hl>;

template <class Par, class M>
using Cubismnc = cubismnc_impl::Cubismnc<Par, M>;
