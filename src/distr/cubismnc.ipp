// Created by Petr Karnakov on 25.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <cassert>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "cubismnc.h"
#include "distr.h"
#include "dump/dumper.h"
#include "util/format.h"
#include "util/histogram.h"

#include "CubismNoCopy/BlockInfo.h"
#include "CubismNoCopy/BlockLab.h"
#include "CubismNoCopy/BlockLabMPI.h"
#include "CubismNoCopy/Grid.h"
#include "CubismNoCopy/GridMPI.h"
#include "CubismNoCopy/HDF5Dumper_MPI.h"
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

  // Convert Cubism BlockInfo to BlockInfoProxy
  static std::vector<BlockInfoProxy> Convert(
      const std::vector<BlockInfo>&, MIdx bs, size_t hl);

  std::vector<size_t> TransferHalos() override;
  std::vector<size_t> TransferHalos(bool inner) override;
  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
  void ReadBuffer(const std::vector<size_t>& bb) override;
  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
  void WriteBuffer(const std::vector<size_t>& bb) override;
  void ReduceSingleRequest(
      const std::vector<std::shared_ptr<RedOp>>& blocks) override;
  void Bcast(const std::vector<size_t>& bb) override;
  void Scatter(const std::vector<size_t>& bb) override;
  void DumpWrite(const std::vector<size_t>& bb) override;

  using P::blocksize_;
  using P::comm_;
  using P::dim;
  using P::extent_;
  using P::frame_;
  using P::halos_;
  using P::isroot_;
  using P::kernelfactory_;
  using P::kernels_;
  using P::nblocks_;
  using P::nprocs_;
  using P::samp_;
  using P::stage_;
  using P::var;

  Grid g_;
  std::vector<std::vector<FieldView>> fviews_; // fields from last communication
  Synch* sync_;
  size_t n_fields_;
  std::map<MIdx, size_t, typename MIdx::LexLess> midx_to_kernel_;
};

// nblocks_ - instance of GFieldViewRaw
template <class nblocks_>
struct StreamHdfScal {
  using B = nblocks_;
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

// nblocks_ - instance of GFieldViewRaw
template <class nblocks_>
struct StreamHdfVect {
  using B = nblocks_;
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
std::vector<BlockInfoProxy> Cubismnc<Par, M>::Convert(
    const std::vector<BlockInfo>& infos, MIdx bs, size_t halos) {
  std::vector<BlockInfoProxy> proxies;
  for (size_t i = 0; i < infos.size(); i++) {
    const BlockInfo& info = infos[i];
    BlockInfoProxy proxy;
    for (int j = 0; j < 3; ++j) {
      proxy.index[j] = info.index[j];
      proxy.origin[j] = info.origin[j];
      proxy.bs[j] = bs[j];
    }
    proxy.h_gridpoint = info.h_gridpoint;
    proxy.hl = halos;
    proxy.maxcomm = 65536;
    proxies.push_back(proxy);
  }
  return proxies;
}

template <class Par, class M>
Cubismnc<Par, M>::Cubismnc(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var)
    : DistrMesh<M>(comm, kf, var)
    , g_(nprocs_[0], nprocs_[1], nprocs_[2], nblocks_[0], nblocks_[1],
         nblocks_[2], extent_, comm) {
  fassert_equal(blocksize_[0], FieldView::bx);
  fassert_equal(blocksize_[1], FieldView::by);
  fassert_equal(blocksize_[2], FieldView::bz);

  int commsize; // size of communicator
  MPI_Comm_size(comm, &commsize);
  if (commsize != nprocs_[0] * nprocs_[1] * nprocs_[2]) {
    throw std::runtime_error(util::Format(
        "Number of MPI tasks ({}) does not match the number of subdomains "
        "px={} py={} pz={}",
        commsize, nprocs_[0], nprocs_[1], nprocs_[2]));
  }

  {
    int rank;
    MPI_Comm_rank(comm, &rank);
    isroot_ = (0 == rank); // XXX: overwrite isroot_
  }

  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Get rid of BlockInfo type
  const std::vector<BlockInfo> infos = g_.getBlocksInfo();
  std::vector<BlockInfoProxy> proxies = Convert(infos, blocksize_, halos_);

  bool islead = true;
  for (size_t i = 0; i < proxies.size(); ++i) {
    auto& proxy = proxies[i];
    const MIdx midx(proxy.index);
    proxy.isroot = (midx == MIdx(0));
    if (proxy.isroot) {
      fassert_equal(i, 0, "Root block must be first");
    }
    proxy.islead = islead;
    const MIdx globalsize = nprocs_ * nblocks_ * blocksize_;
    for (int j = 0; j < 3; ++j) {
      proxy.gs[j] = globalsize[j];
    }
    islead = false;
    std::cerr << util::Format("q=({} {})", midx, i) << " ";
    midx_to_kernel_[midx] = i;
  }
  comm_ = g_.getCartComm(); // XXX: overwrite comm_
  this->MakeKernels(proxies);
}

template <class Par, class M>
Cubismnc<Par, M>::~Cubismnc() {
  for (const auto& sm : g_.getSynchronizerMPI()) {
    samp_.Append(sm.second->getSampler());
  }
}

template <class Par, class M>
auto Cubismnc<Par, M>::TransferHalos() -> std::vector<size_t> {
  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] not needed
  samp_.SeedSample();
  MPI_Barrier(comm_);
  samp_.CollectSample("MPI_Barrier");

  const std::vector<BlockInfo> infos = g_.getBlocksInfo(); // all blocks
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
      assert(n_fields == fields.size());
      for (size_t fi = 0; fi < n_fields; ++fi) {
        auto& req = reqs[fi];
        fviews[fi].push_back(
            FieldView(req->GetBasePtr(), infos[i].index, req->GetSize()));
      }
    }

    // 1. Exchange halos in buffer mesh.
    // stencil type
    const bool is_tensorial = true;
    const int nhalo_start[3] = {halos_, halos_, halos_};
    const int nhalo_end[3] = {halos_, halos_, halos_};

    // schedule asynchronous communication
    Synch& s = g_.sync(
        fviews, nhalo_start, nhalo_end, is_tensorial,
        var.Int["mpi_compress_msg"], var.Int["histogram"]);

    // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] not needed
    samp_.SeedSample();
    MPI_Barrier(comm_);
    samp_.CollectSample("MPI_Barrier");

    // Get all blocks synchronized
    const std::vector<BlockInfo> avail = s.avail_halo();
    fassert_equal(avail.size(), infos.size());

    // 2. Load exchanged halos into the local fields
    // FIXME: [fabianw@mavt.ethz.ch; 2019-12-07] This parallel structure is
    // likely less efficient for single threaded execution
#pragma omp parallel
    {
      Lab l;
      l.prepare(g_, s);
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
  const std::vector<BlockInfo> infos = g_.getBlocksInfo(); // all blocks
  fassert_equal(infos.size(), kernels_.size());

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
      const int nhalo_start[3] = {halos_, halos_, halos_};
      const int nhalo_end[3] = {halos_, halos_, halos_};

      // schedule asynchronous communication
      sync_ = &g_.sync(
          fviews_, nhalo_start, nhalo_end, is_tensorial,
          var.Int["mpi_compress_msg"], var.Int["histogram"]);
      const std::vector<BlockInfo> avail = sync_->avail_inner();

      {
        std::cerr << "*1***\n";
        for (auto& info : avail) {
          std::cerr << " " << MIdx(info.index);
        }
        std::cerr << "\n";
      }
      {
        std::cerr << "*midx_to_kernel***\n";
        for (auto q : midx_to_kernel_) {
          std::cerr << util::Format("({} {})", q.first, q.second) << " ";
        }
        std::cerr << "\n";
      }

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
        l.prepare(g_, *sync_);
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
        l.prepare(g_, *sync_);
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

  std::cerr << util::Format(
                   "stage={} bb={} isroot_={} n_fields={} inner={}", stage_, bb,
                   isroot_, n_fields_, inner)
            << std::endl;
  return bb;
}

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
template <class Par, class M>
void Cubismnc<Par, M>::ReadBuffer(const std::vector<size_t>&) {}

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
template <class Par, class M>
void Cubismnc<Par, M>::WriteBuffer(const std::vector<size_t>&) {}

template <class Par, class M>
void Cubismnc<Par, M>::Bcast(const std::vector<size_t>& bb) {
  using OpConcat = typename M::OpCat;
  auto& vf = kernels_.front()->GetMesh().GetBcast(); // pointers to broadcast

  // Check size is the same for all kernels
  for (auto b : bb) {
    auto& v = kernels_[b]->GetMesh().GetBcast(); // pointers to broadcast
    fassert_equal(v.size(), vf.size());
  }

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpConcat* o = dynamic_cast<OpConcat*>(vf[i].get())) {
      std::vector<char> r = o->Neutral(); // buffer

      if (isroot_) {
        // read from root block
        for (auto b : bb) {
          auto& m = kernels_[b]->GetMesh();
          if (m.IsRoot()) {
            auto& v = m.GetBcast();
            OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
            ob->Append(r);
          }
        }
      }

      int s = r.size(); // size

      // broadcast size
      MPI_Bcast(&s, 1, MPI_INT, 0, comm_);

      // resize
      r.resize(s);

      // broadcast data
      samp_.SeedSample();
      MPI_Bcast(r.data(), r.size(), MPI_CHAR, 0, comm_);
      samp_.CollectSample("MPI_Bcast");

      // write to all blocks
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        auto& v = m.GetBcast();
        OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
        ob->Set(r);
      }
    } else {
      throw std::runtime_error("Bcast: Unknown M::Op instance");
    }
  }

  // Clear bcast requests
  for (auto b : bb) {
    auto& k = *kernels_[b];
    auto& m = k.GetMesh();
    m.ClearBcast();
  }
}
template <class Par, class M>
void Cubismnc<Par, M>::Scatter(const std::vector<size_t>& bb) {
  auto& vreq0 =
      kernels_.front()->GetMesh().GetScatter(); // requests on first block

  // Check size is the same for all blocks
  for (auto b : bb) {
    auto& vreq = kernels_[b]->GetMesh().GetScatter();
    fassert_equal(vreq.size(), vreq0.size());
  }

  for (size_t q = 0; q < vreq0.size(); ++q) {
    int recvcount;
    int sizes_recvcount;
    std::vector<Scal> rbuf;
    std::vector<int> sizes_rbuf;

    MPI_Datatype mscal = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    if (isroot_) {
      int commsize; // size of communicator
      MPI_Comm_size(comm_, &commsize);
      std::vector<Scal> buf;
      std::vector<int> dis(commsize, 0);
      std::vector<int> cnt(commsize, 0);
      std::vector<int> sizes_buf;
      std::vector<int> sizes_dis(commsize, 0);
      std::vector<int> sizes_cnt(commsize, 0);
      // find root block
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        if (m.IsRoot()) {
          auto& req = m.GetScatter()[q];
          size_t i = 0;
          // concatenate data for all blocks in buf
          for (int rank = 0; rank < commsize; ++rank) {
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
      samp_.SeedSample();
      MPI_Scatterv(
          buf.data(), cnt.data(), dis.data(), mscal, rbuf.data(), recvcount,
          mscal, 0, comm_);
      samp_.CollectSample("MPI_Scatterv");
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
      samp_.SeedSample();
      MPI_Scatterv(
          nullptr, nullptr, nullptr, mscal, rbuf.data(), recvcount, mscal, 0,
          comm_);
      samp_.CollectSample("MPI_Scatterv");
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
void Cubismnc<Par, M>::ReduceSingleRequest(
    const std::vector<std::shared_ptr<RedOp>>& blocks) {
  using OpScal = typename M::OpS;
  using OpScalInt = typename M::OpSI;
  using OpConcat = typename M::OpCat;

  auto* firstbase = blocks.front().get();

  if (auto* first = dynamic_cast<OpScal*>(firstbase)) {
    auto buf = first->Neutral(); // result

    // Reduce over blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase.get());
      other->Append(buf);
    }

    MPI_Op mpiop;
    if (dynamic_cast<typename M::OpSum*>(first)) {
      mpiop = MPI_SUM;
    } else if (dynamic_cast<typename M::OpProd*>(first)) {
      mpiop = MPI_PROD;
    } else if (dynamic_cast<typename M::OpMax*>(first)) {
      mpiop = MPI_MAX;
    } else if (dynamic_cast<typename M::OpMin*>(first)) {
      mpiop = MPI_MIN;
    } else {
      fassert(false, "Unknown reduction");
    }
    MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    // Reduce over ranks
    MPI_Allreduce(MPI_IN_PLACE, &buf, 1, mt, mpiop, comm_);

    // Write results to all blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase.get());
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpScalInt*>(firstbase)) {
    auto buf = first->Neutral(); // result

    // Reduce over blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase.get());
      other->Append(buf);
    }

    MPI_Op mpiop;
    if (dynamic_cast<typename M::OpMinloc*>(first)) {
      mpiop = MPI_MINLOC;
    } else if (dynamic_cast<typename M::OpMaxloc*>(first)) {
      mpiop = MPI_MAXLOC;
    } else {
      fassert(false, "Unknown reduction");
    }

    MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE_INT : MPI_FLOAT_INT);

    // Reduce over all ranks
    MPI_Allreduce(MPI_IN_PLACE, &buf, 1, mt, mpiop, comm_);

    // Write results to all blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase.get());
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpConcat*>(firstbase)) {
    auto buf = first->Neutral();

    // Reduce over local blocks
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpConcat*>(otherbase.get());
      other->Append(buf);
    }

    int bufsize = buf.size();

    if (isroot_) {
      int commsize;
      MPI_Comm_size(comm_, &commsize);

      std::vector<int> sizes(commsize);

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
    if (dumpformat == "default") {
      dumpformat = "hdf";
    }

    if (dumpformat == "hdf") {
      // Create FieldView's for dump
      const size_t n_fields = mfirst.GetDump().size();
      std::vector<BlockInfo> infos = g_.getBlocksInfo(); // all blocks
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
              blocks, aos_idx, g_, frame_, frame_, name, ".", Vect(0),
              mfirst.GetCellSize(), true);
        } else if (3 == blocks[0].n_comp) {
          StreamHdfVect<FieldView>::NAME = reqs[fi].second;
          DumpHDF5_MPI<StreamHdfVect<FieldView>>(
              blocks, aos_idx, g_, frame_, frame_, name, ".", Vect(0),
              mfirst.GetCellSize(), true);
        } else {
          throw std::runtime_error("DumpWrite(): Support only size 1 and 3");
        }
      }
      if (isroot_) {
        std::cout << "Dump " << frame_ << ": format=" << dumpformat
                  << std::endl;
      }
      ++frame_;
    } else {
      P::DumpWrite(bb);
    }
  }
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
