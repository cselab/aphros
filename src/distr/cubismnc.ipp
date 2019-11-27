// vim: expandtab:smarttab:sw=2:ts=2
#pragma once

#include <cassert>
#include <memory>
#include <map>
#include <stdexcept>
#include <mpi.h>

#include "distr.h"
#include "cubismnc.h"
#include "dump/dumper.h"

#include "CubismNoCopy/BlockLab.h"
#include "CubismNoCopy/BlockLabMPI.h"
#include "CubismNoCopy/Grid.h"
#include "CubismNoCopy/GridMPI.h"
#include "CubismNoCopy/HDF5Dumper_MPI.h"
#include "CubismNoCopy/BlockInfo.h"
#include "CubismNoCopy/StencilInfo.h"

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
    GFieldViewRaw(Elem *const base, const int idx[3], const int n_components=1)
        : data(base + n_components * (halo + halo * stridex + halo * stridex
                                      * stridey)),
          data_halo(base), index{idx[0], idx[1], idx[2]}, n_comp(n_components)
    {
      assert(data_halo != data); // halo = 0 will not work in this model
    }

    // SoA field
    Elem *const data;      // block data
    Elem *const data_halo; // block data including halos

    const int index[3]; // corresponding block index

    // Number of components carried by field.  Parameter is not passed as
    // template parameter such that scalar and vector (field) views can be
    // combined in the same field compound (to maximize MPI communication buffer
    // size).
    const int n_comp;

    // block access emulator
    inline Elem &operator()(unsigned int x, unsigned int y = 0, unsigned int z = 0)
    {
        assert(data != nullptr);
        assert(0 <= x && x < (int)bx);
        assert(0 <= y && y < (int)by);
        assert(0 <= z && z < (int)bz);
        return data[n_comp * (x + stridex * (y + stridey * z))];
    }

    inline Elem &LinAccess(unsigned int x)
    {
        assert(data_halo != nullptr);
        assert(x < n_comp * stridex * stridey * stridez);

        return data_halo[n_comp * x];
    }
};

// TView - instance of GFieldViewRaw
template <class TView>
class LabPer : public BlockLab<TView>
{
    using FieldView = TView;
    using ElementTypeBlock = typename FieldView::Elem;

public:
    virtual inline std::string name() const { return "LabPer"; }
    bool is_xperiodic() { return true; }
    bool is_yperiodic() { return true; }
    bool is_zperiodic() { return true; }
};

// TView - instance of GFieldViewRaw
template <class TView>
using GLab = BlockLabMPI<LabPer<TView>>;

// TView - instance of GFieldViewRaw
template <class TView>
using GGrid = GridMPI<Grid<TView>>;

// Par - instance of GPar
// KF - instance of KernelMeshFactory
template <class Par, class KF>
class Cubismnc : public DistrMesh<KF> {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  Cubismnc(MPI_Comm comm, KF& kf, Vars& var);
  typename M::BlockCells GetGlobalBlock() const override;
  typename M::IndexCells GetGlobalIndex() const override;
  FieldCell<Scal> GetGlobalField(size_t i) override;

 private:
  using FieldView = GFieldViewRaw<Par>;
  using Lab = GLab<FieldView>;
  using Grid = GGrid<FieldView>;
  using Elem = typename FieldView::Elem;
  using Synch = typename Grid::Synch;

  using P = DistrMesh<KF>; // parent
  using MIdx = typename M::MIdx;

  using P::dim;
  using P::mk;
  using P::kf_;
  using P::var;
  using P::bs_;
  using P::hl_;
  using P::p_;
  using P::b_;
  using P::stage_;
  using P::isroot_;
  using P::comm_;
  using P::ext_;
  using P::frame_;

  Grid g_;
// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] The map is not really needed
  struct S { // cubism [s]tate
    std::map<MIdx, BlockInfo, typename MIdx::LexLess> mb;
  };
  S s_;
  std::vector<std::vector<FieldView>> fviews_; // fields from last communication
  Synch* sync_;
  size_t n_fields_;

  // Convert Cubism BlockInfo to MyBlockInfo
  static std::vector<MyBlockInfo> Convert(
      const std::vector<BlockInfo>&, MIdx bs, size_t hl);

  std::vector<MIdx> GetBlocks() override;
  std::vector<MIdx> GetBlocks(bool inner) override;
// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
  void ReadBuffer(const std::vector<MIdx>& bb) override;
// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void Reduce(const std::vector<MIdx>& bb) override;
  void Bcast(const std::vector<MIdx>& bb) override;
  void Scatter(const std::vector<MIdx>& bb) override;
  void DumpWrite(const std::vector<MIdx>& bb) override;
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

  StreamHdfScal(B &b, const int idx) : b(b), aos(idx)
  {
      assert(1 == b.n_comp);
      assert(0 >= aos);
  }

  // write
  void operate(const int ix, const int iy, const int iz, Scal* out) const
  {
    out[0] = (&b(ix, iy, iz))[aos];
  }

  // read
  void operate(const Scal* out, const int ix, const int iy, const int iz)
  {
    (&b(ix, iy, iz))[aos] = out[0];
  }

  static const char * getAttributeName() { return "Scalar"; }
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

  StreamHdfVect(B &b, const int idx) : b(b), aos(idx)
  {
      assert(3 == b.n_comp);
      assert(0 == aos);
  }

  // write
  void operate(const int ix, const int iy, const int iz, Scal* out) const
  {
    const Scal *v = &((const Scal *)&b(ix, iy, iz))[aos];
    out[0] = v[0];
    out[1] = v[1];
    out[2] = v[2];
  }

  // read
  void operate(const Scal* out, const int ix, const int iy, const int iz)
  {
    Scal *v = &((Scal *)&b(ix, iy, iz))[aos];
    v[0] = out[0];
    v[1] = out[1];
    v[2] = out[2];
  }

  static const char * getAttributeName() { return "Vector"; }
};

template <class Par, class KF>
std::vector<MyBlockInfo> Cubismnc<Par, KF>::Convert(
    const std::vector<BlockInfo>& cc, MIdx bs, size_t hl) {
  std::vector<MyBlockInfo> bb;
  for(size_t i = 0; i < cc.size(); i++) {
    const BlockInfo& c = cc[i];
    MyBlockInfo b;
    for (int j = 0; j < 3; ++j) {
      b.index[j] = c.index[j];
      b.origin[j] = c.origin[j];
      b.bs[j] = bs[j];
    }
    b.h_gridpoint = c.h_gridpoint;
    b.hl = hl;
    b.maxcomm = 65536;
    bb.push_back(b);
  }
  return bb;
}


template <class Par, class KF>
Cubismnc<Par, KF>::Cubismnc(MPI_Comm comm, KF& kf, Vars& var)
  : DistrMesh<KF>(comm, kf, var)
  , g_(p_[0], p_[1], p_[2], b_[0], b_[1], b_[2], ext_, comm)
{
  assert(bs_[0] == FieldView::bx &&
         bs_[1] == FieldView::by &&
         bs_[2] == FieldView::bz);

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);  // XXX: overwrite isroot_

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Get rid of BlockInfo type
  std::vector<BlockInfo> cc = g_.getBlocksInfo();
  std::vector<MyBlockInfo> ee = Convert(cc, bs_, hl_);

  bool islead = true;
  for (auto& e : ee) {
    MIdx d(e.index);
    e.isroot = (d == MIdx(0));
    e.islead = islead;
    MIdx gs = p_ * b_ * bs_; // global size
    for (int j = 0; j < 3; ++j) {
      e.gs[j] = gs[j];
    }
    islead = false;
  }

  this->MakeKernels(ee);

  comm_ = g_.getCartComm(); // XXX: overwrite comm_
}

template <class Par, class KF>
auto Cubismnc<Par, KF>::GetBlocks() -> std::vector<MIdx> {
  // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] not needed
  MPI_Barrier(comm_);

  // Get all blocks
  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // all blocks
  auto& m = mk.at(MIdx(cc[0].index))->GetMesh();
  std::vector<BlockInfo> aa;
  // Perform communication if necessary
  if (m.GetComm().size() > 0) {
    // 0. Construct vector of fields to communicate for all blocks on this rank
    std::vector<std::vector<FieldView>> fviews;
    const size_t n_fields = m.GetComm().size();
    fviews.resize(n_fields);
    for (const auto &b : cc) {
      auto &bm = mk.at(MIdx(b.index))->GetMesh();
      auto &bf = bm.GetComm();
      assert(n_fields == bf.size());
      for (size_t fi = 0; fi < n_fields; ++fi) {
        auto &o = bf[fi];
        fviews[fi].push_back(FieldView(o->GetBasePtr(), b.index, o->GetSize()));
      }
    }

    // 1. Exchange halos in buffer mesh.
    // stencil type
    const bool is_tensorial = true;
    const int nhalo_start[3] = {hl_, hl_, hl_};
    const int nhalo_end[3] = {hl_, hl_, hl_};

    // schedule asynchronous communication
    Synch& s = g_.sync(fviews, nhalo_start, nhalo_end, is_tensorial);

    // FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] not needed
    MPI_Barrier(comm_);

    // Get all blocks synchronized
    aa = s.avail(cc.size());

    // 2. Load exchanged halos into the local fields
    Lab l;
    l.prepare(g_, s);
    for (auto& field : fviews) {
      for (auto& block : field) {
        l.load(block, field);
      }
    }
  } else {
    aa = cc;
  }

  if (aa.size() != cc.size()) {
    std::cerr
        << "aa.size()=" << aa.size()
        << " != "
        << "cc.size()=" << cc.size()
        << std::endl;
    assert(false);
  }

  // Create vector of indices and save block info to map
  std::vector<MIdx> bb;
  s_.mb.clear();
  for (auto a : aa) {
    MIdx b(a.index);
    s_.mb.emplace(b, a);
    bb.push_back(b);
  }

  return bb;
}

template <class Par, class KF>
auto Cubismnc<Par, KF>::GetBlocks(bool inner) -> std::vector<MIdx> {
  if (inner) {
    std::vector<BlockInfo> cc = g_.getBlocksInfo(); // all blocks
    n_fields_ = mk.at(MIdx(cc[0].index))->GetMesh().GetComm().size();
    std::vector<MIdx> bb;
    // Perform communication if necessary
    if (n_fields_ > 0) {
      // 0. Construct vector of fields to communicate for all blocks on this rank
      fviews_.clear();
      fviews_.resize(n_fields_);
      for (const auto &b : cc) {
        auto &bm = mk.at(MIdx(b.index))->GetMesh();
        auto &bf = bm.GetComm();
        assert(n_fields_ == bf.size());
        for (size_t fi = 0; fi < n_fields_; ++fi) {
          auto &o = bf[fi];
          fviews_[fi].push_back(FieldView(o->GetBasePtr(), b.index, o->GetSize()));
        }
      }

      // 1. Exchange halos in buffer mesh.
      // stencil type
      const bool is_tensorial = true;
      const int nhalo_start[3] = {hl_, hl_, hl_};
      const int nhalo_end[3] = {hl_, hl_, hl_};

      // schedule asynchronous communication
      sync_ = &g_.sync(fviews_, nhalo_start, nhalo_end, is_tensorial);
      std::vector<BlockInfo> aa = sync_->avail_inner();

      // Create vector of indices and save block info to map
      s_.mb.clear();
      for (auto a : aa) {
        MIdx b(a.index);
        s_.mb.emplace(b, a);
        bb.push_back(b);
      }

      std::set<MIdx> set(bb.begin(), bb.end());

      // 2. Load exchanged halos into the local fields
      Lab l;
      l.prepare(g_, *sync_);
      for (auto& field : fviews_) {
        for (auto& block : field) {
          if (set.count(MIdx(block.index))) {
            l.load(block, field);
          }
        }
      }
    } else {
      s_.mb.clear();
      for (auto a : cc) {
        MIdx b(a.index);
        s_.mb.emplace(b, a);
        bb.push_back(b);
      }
    }
    return bb;
  } else {
    std::vector<MIdx> bb;
    // Perform communication if necessary
    if (n_fields_ > 0) {
      std::vector<BlockInfo> aa = sync_->avail_halo();

      for (auto a : aa) {
        MIdx b(a.index);
        s_.mb.emplace(b, a);
        bb.push_back(b);
      }

      std::set<MIdx> set(bb.begin(), bb.end());

      // Load exchanged halos into the local fields
      Lab l;
      l.prepare(g_, *sync_);
      for (auto& field : fviews_) {
        for (auto& block : field) {
          if (set.count(MIdx(block.index))) {
            l.load(block, field);
          }
        }
      }
    }
    return bb;
  }
}

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
template <class Par, class KF>
void Cubismnc<Par, KF>::ReadBuffer(const std::vector<MIdx>&) {}

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
template <class Par, class KF>
void Cubismnc<Par, KF>::WriteBuffer(const std::vector<MIdx>&) {}

template <class Par, class KF>
void Cubismnc<Par, KF>::Bcast(const std::vector<MIdx>& bb) {
  using OpCat = typename M::OpCat;
  auto& vf = mk.at(bb[0])->GetMesh().GetBcast();  // pointers to broadcast

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& v = mk.at(b)->GetMesh().GetBcast();  // pointers to broadcast
    if (v.size() != vf.size()) {
      throw std::runtime_error("Bcast: v.size() != vf.size()");
    }
  }

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpCat* o = dynamic_cast<OpCat*>(vf[i].get())) {
      std::vector<char> r = o->Neut(); // buffer

      if (isroot_) {
        // read from root block
        for (auto& b : bb) {
          auto& m = mk.at(b)->GetMesh();
          if (m.IsRoot()) {
            auto& v = m.GetBcast();
            OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
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
      MPI_Bcast(r.data(), r.size(), MPI_CHAR, 0, comm_);

      // write to all blocks
      for (auto& b : bb) {
        auto& m = mk.at(b)->GetMesh();
        auto& v = m.GetBcast();
        OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
        ob->Set(r);
      }
    } else {
      throw std::runtime_error("Bcast: Unknown M::Op instance");
    }
  }

  // Clear bcast requests
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    auto& m = k.GetMesh();
    m.ClearBcast();
  }
}
template <class Par, class KF>
void Cubismnc<Par, KF>::Scatter(const std::vector<MIdx>& bb) {
  auto& vreq0 = mk.at(bb[0])->GetMesh().GetScatter(); // requests on first block

  // Check size is the same for all blocks
  for (auto& b : bb) {
    auto& vreq = mk.at(b)->GetMesh().GetScatter();
    if (vreq.size() != vreq0.size()) {
      throw std::runtime_error("Scatter: vreq.size() != vreq0.size()");
    }
  }

  for (size_t q = 0; q < vreq0.size(); ++q) {
    int recvcount;
    int sizes_recvcount;
    std::vector<Scal> rbuf;
    std::vector<int> sizes_rbuf;

    MPI_Datatype mscal = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    if (isroot_) {
      int sc; // size of communicator
      MPI_Comm_size(comm_, &sc);
      std::vector<Scal> buf;
      std::vector<int> dis(sc, 0);
      std::vector<int> cnt(sc, 0);
      std::vector<int> sizes_buf;
      std::vector<int> sizes_dis(sc, 0);
      std::vector<int> sizes_cnt(sc, 0);
      // find root block
      for (auto& b : bb) {
        auto& m = mk.at(b)->GetMesh();
        if (m.IsRoot()) {
          auto& req = m.GetScatter()[q];
          size_t i = 0;
          // concatenate data for all blocks in buf
          for (int rank = 0; rank < sc; ++rank) {
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
      MPI_Scatter(cnt.data(), 1, MPI_INT,
                  &recvcount, 1, MPI_INT, 0, comm_);
      rbuf.resize(recvcount);
      // data
      MPI_Scatterv(buf.data(), cnt.data(), dis.data(), mscal,
                   rbuf.data(), recvcount, mscal, 0, comm_);
      // sizes recvcount
      MPI_Scatter(sizes_cnt.data(), 1, MPI_INT,
                  &sizes_recvcount, 1, MPI_INT, 0, comm_);
      sizes_rbuf.resize(sizes_recvcount);
      // sizes
      MPI_Scatterv(
          sizes_buf.data(), sizes_cnt.data(), sizes_dis.data(), MPI_INT,
          sizes_rbuf.data(), sizes_recvcount, MPI_INT, 0, comm_);
    } else {
      // data recvcount
      MPI_Scatter(nullptr, 0, MPI_INT,
                  &recvcount, 1, MPI_INT, 0, comm_);
      rbuf.resize(recvcount);
      // data
      MPI_Scatterv(nullptr, nullptr, nullptr, mscal,
                   rbuf.data(), recvcount, mscal, 0, comm_);
      // sizes recvcount
      MPI_Scatter(nullptr, 0, MPI_INT,
                  &sizes_recvcount, 1, MPI_INT, 0, comm_);
      sizes_rbuf.resize(sizes_recvcount);
      // sizes
      MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT,
                   sizes_rbuf.data(), sizes_recvcount, MPI_INT, 0, comm_);
    }

    // write to blocks on current rank
    size_t off = 0;
    for (size_t k = 0; k < bb.size(); ++k) {
      auto& v = *mk.at(bb[k])->GetMesh().GetScatter()[q].second;
      v = std::vector<Scal>(rbuf.data() + off,
                            rbuf.data() + off + sizes_rbuf[k]);
      off += sizes_rbuf[k];
    }
  }

  // Clear requests
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearScatter();
  }
}
template <class Par, class KF>
void Cubismnc<Par, KF>::Reduce(const std::vector<MIdx>& bb) {
  using OpS = typename M::OpS;
  using OpSI = typename M::OpSI;
  using OpCat = typename M::OpCat;
  auto& vf = mk.at(bb[0])->GetMesh().GetReduce();  // pointers to reduce

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& v = mk.at(b)->GetMesh().GetReduce();  // pointers to reduce
    if (v.size() != vf.size()) {
      throw std::runtime_error("Reduce: v.size() != vf.size()");
    }
  }

  // TODO: Check operation is the same for all kernels
  // TODO: avoid code duplication

  for (size_t i = 0; i < vf.size(); ++i) {
    if (OpS* o = dynamic_cast<OpS*>(vf[i].get())) {
      auto r = o->Neut(); // result

      // Reduce over all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce();
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Append(r);
      }

      MPI_Op mo;
      if (dynamic_cast<typename M::OpSum*>(o)) {
        mo = MPI_SUM;
      } else if (dynamic_cast<typename M::OpProd*>(o)) {
        mo = MPI_PROD;
      } else if (dynamic_cast<typename M::OpMax*>(o)) {
        mo = MPI_MAX;
      } else if (dynamic_cast<typename M::OpMin*>(o)) {
        mo = MPI_MIN;
      } else {
        throw std::runtime_error("Reduce(): Can't find MPI_Op");
      }
      MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

      // Reduce over all ranks
      MPI_Allreduce(MPI_IN_PLACE, &r, 1, mt, mo, comm_);

      // Write results to all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce();
        OpS* ob = dynamic_cast<OpS*>(v[i].get());
        ob->Set(r);
      }
    } else if (OpSI* o = dynamic_cast<OpSI*>(vf[i].get())) {
      auto r = o->Neut(); // result

      // Reduce over all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce();
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Append(r);
      }

      MPI_Op mo;
      if (dynamic_cast<typename M::OpMinloc*>(o)) {
        mo = MPI_MINLOC;
      } else if (dynamic_cast<typename M::OpMaxloc*>(o)) {
        mo = MPI_MAXLOC;
      } else {
        throw std::runtime_error("Reduce(): Can't find MPI_Op");
      }

      MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE_INT : MPI_FLOAT_INT);

      // Reduce over all ranks
      MPI_Allreduce(MPI_IN_PLACE, &r, 1, mt, mo, comm_);

      // Write results to all blocks on current rank
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce();
        OpSI* ob = dynamic_cast<OpSI*>(v[i].get());
        ob->Set(r);
      }
    } else if (OpCat* o = dynamic_cast<OpCat*>(vf[i].get())) {
      std::vector<char> r = o->Neut(); // result local

      // Reduce over all local blocks
      for (auto& b : bb) {
        auto& v = mk.at(b)->GetMesh().GetReduce();
        OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
        ob->Append(r);
      }

      int s = r.size(); // size local

      if (isroot_) {
        int sc; // size of communicator
        MPI_Comm_size(comm_, &sc);

        std::vector<int> ss(sc); // size of r on all ranks

        // Gather ss
        MPI_Gather(&s, 1, MPI_INT,
                   ss.data(), 1, MPI_INT,
                   0, comm_);


        int sa = 0; // size all
        std::vector<int> oo = {0}; // offsets
        for (auto& q : ss) {
          sa += q;
          oo.push_back(oo.back() + q);
        }
        oo.pop_back();
        assert(ss.size() == oo.size());

        std::vector<char> ra(sa); // result all

        // Gather ra
        MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                    ra.data(), ss.data(), oo.data(), MPI_CHAR,
                    0, comm_);

        // Write results to root block
        size_t cnt = 0;
        for (auto& b : bb) {
          auto& m = mk.at(b)->GetMesh();
          if (m.IsRoot()) {
            auto& v = m.GetReduce();
            OpCat* ob = dynamic_cast<OpCat*>(v[i].get());
            ob->Set(ra);
            ++cnt;
          }
        }
        assert(cnt == 1);
      } else {
        // Send s to root
        MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm_);

        // Send r to root
        MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                    nullptr, nullptr, nullptr, MPI_CHAR,
                    0, comm_);
      }

    } else {
      throw std::runtime_error("Reduce: Unknown M::Op implementation");
    }
  }

  // Clear reduce requests
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    auto& m = k.GetMesh();
    m.ClearReduce();
  }
}

template <class Par, class KF>
void Cubismnc<Par, KF>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = var.String["dumpformat"];
    if (df == "default") {
      df = "hdf";
    }

    if (df == "hdf") {
      // Create FieldView's for dump
      const size_t n_fields = m.GetDump().size();
      std::vector<BlockInfo> cc = g_.getBlocksInfo(); // all blocks
      std::vector<std::vector<FieldView>> fviews(n_fields);
      for (const auto &b : cc) {
        auto &bm = mk.at(MIdx(b.index))->GetMesh();
        auto &bf = bm.GetDump();
        assert(n_fields == bf.size());
        for (size_t fi = 0; fi < n_fields; ++fi) {
          auto &o = bf[fi].first;
          fviews[fi].push_back(
            FieldView(o->GetBasePtr(), b.index, o->GetStride()));
        }
      }

      // Write dump
      auto &bf = m.GetDump();
      assert(bf.size() == fviews.size());
      for (size_t fi = 0; fi < fviews.size(); ++fi) {
        auto fn = GetDumpName(bf[fi].second, "", frame_);
        const int aos_idx = (bf[fi].first)->GetIndex();
        auto &blocks = fviews[fi];
        if (0 <= aos_idx) {
          StreamHdfScal<FieldView>::NAME = bf[fi].second;
          DumpHDF5_MPI<StreamHdfScal<FieldView>>(blocks, aos_idx, g_, frame_,
            frame_, fn, ".", Vect(0), m.GetCellSize(), true);
        } else if (3 == blocks[0].n_comp) {
          StreamHdfVect<FieldView>::NAME = bf[fi].second;
          DumpHDF5_MPI<StreamHdfVect<FieldView>>(blocks, aos_idx, g_, frame_,
            frame_, fn, ".", Vect(0), m.GetCellSize(), true);
        } else {
          throw std::runtime_error("DumpWrite(): Support only size 1 and 3");
        }
      }
      if (isroot_) {
        std::cerr << "Dump " << frame_ << ": format=" << df << std::endl;
      }
      ++frame_;
    } else {
      P::DumpWrite(bb);
    }
  }
}

template <class Par, class KF>
auto Cubismnc<Par, KF>::GetGlobalBlock() const -> typename M::BlockCells {
  return typename M::BlockCells(p_ * b_ * bs_);
}

template <class Par, class KF>
auto Cubismnc<Par, KF>::GetGlobalIndex() const -> typename M::IndexCells {
  // TODO: revise, relies on lightweight mesh,
  //       may cause allocation of data fields otherwise
  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx s = GetGlobalBlock().GetSize();
  auto m = InitUniformMesh<M>(dom, MIdx(0), s, 0, true, true, s, 0);
  return m.GetIndexCells();
}

template <class Par, class KF>
auto Cubismnc<Par, KF>::GetGlobalField(size_t e) -> FieldCell<Scal> {
  using BC = typename M::BlockCells;
  using PF = const typename M::Co*;
  auto gbc = GetGlobalIndex();
  // collective, does actual communication
  auto bb = GetBlocks();
  std::vector<Scal> v(bs_.prod()); // tmp
  MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);
  BC bc(bs_); // cells of one block
  GBlock<size_t, dim> bq(p_ * b_);  // indices of block
  GIndex<size_t, dim> ndq(p_ * b_);  // flat
  auto GetField = [e](const M& m) -> PF {
    int k = static_cast<int>(e);
    for (auto& c : m.GetComm()) { // must be first
      k -= c->GetSize();
      if (k < 0) {
        return c.get();
      }
    }
    for (auto& co : m.GetDump()) { // followed by this
      k -= co.first->GetSize();
      if (k < 0) {
        return co.first.get();
      }
    }
    return nullptr;
  };
  if (isroot_) {
    FieldCell<Scal> gfc(gbc); // result
    // Copy from blocks on root
    for (auto& b : bb) {
      // block mesh
      auto& m = mk.at(b)->GetMesh();
      // index cells
      auto& mbc = m.GetIndexCells();
      // get corner of inner cells block
      MIdx wb = m.GetInBlockCells().GetBegin();
      // get field fc associated to index e
      const PF bf = GetField(m);
      if (const auto fcptr = dynamic_cast<const typename M::CoFcs*>(bf)) {
        const FieldCell<Scal> &fc = *(fcptr->f);
        // copy from inner cells to global field
        for (auto w : bc) {
          gfc[gbc.GetIdx(wb + w)] = fc[mbc.GetIdx(wb + w)];
        }
      } else if (const auto fcptr = dynamic_cast<const typename M::CoFcv*>(bf)) {
        const FieldCell<Vect> &fc = *(fcptr->f);
        // copy from inner cells to global field
        for (auto w : bc) {
          gfc[gbc.GetIdx(wb + w)] = fc[mbc.GetIdx(wb + w)][fcptr->d];
        }
      } else {
        throw std::runtime_error("GetGlobalField: resolving field pointer failed");
      }
    }
    // recv from other ranks
    for (auto b : bq) {
      if (!s_.mb.count(b)) { // not local block
        MPI_Status st;
        MPI_Recv(v.data(), v.size(), mt, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, comm_, &st);

        size_t i = 0;
        MIdx wb = ndq.GetMIdx(size_t(st.MPI_TAG)) * bs_;
        for (auto w : bc) {
          gfc[gbc.GetIdx(wb + w)] = v[i++];
        }
      }
    }

    MPI_Barrier(comm_);
    return gfc;
  } else {
    // send to root
    for (auto b : bb) {
      // block mesh
      auto& m = mk.at(b)->GetMesh();
      // block cells
      auto& mbc = m.GetIndexCells();
      // get corner of inner cells block
      MIdx wb = m.GetInBlockCells().GetBegin();
      // get field fc associated to index e
      const PF bf = GetField(m);
      if (const auto fcptr = dynamic_cast<const typename M::CoFcs*>(bf)) {
        const FieldCell<Scal> &fc = *(fcptr->f);
        size_t i = 0;
        // copy from inner cells to global field
        for (auto w : bc) {
          v[i++] = fc[mbc.GetIdx(wb + w)];
        }
      } else if (const auto fcptr = dynamic_cast<const typename M::CoFcv*>(bf)) {
        const FieldCell<Vect> &fc = *(fcptr->f);
        size_t i = 0;
        // copy from inner cells to global field
        for (auto w : bc) {
          v[i++] = fc[mbc.GetIdx(wb + w)][fcptr->d];
        }
      } else {
        throw std::runtime_error("GetGlobalField: failed dynamic cast");
      }
      // XXX: assume same order of Recv on root
      MPI_Send(v.data(), v.size(), mt, 0, ndq.GetIdx(b), comm_);
    }
    MPI_Barrier(comm_);
    return FieldCell<Scal>();
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

template <class Par, class KF>
using Cubismnc = cubismnc_impl::Cubismnc<Par, KF>;
