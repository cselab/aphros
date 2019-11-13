// vim: expandtab:smarttab:sw=2:ts=2
#pragma once

#include <cassert>
#include <memory>
#include <limits>
#include <map>
#include <stdexcept>
#include <mpi.h>

#include "distr.h"
#include "cubism.h"
#include "dump/dumper.h"

#include "Cubism/afree/BlockLab.h"
#include "Cubism/afree/BlockLabMPI.h"
#include "Cubism/afree/Grid.h"
#include "Cubism/afree/GridMPI.h"
#include "Cubism/afree/HDF5Dumper_MPI.h"
#include "Cubism/BlockInfo.h"
#include "Cubism/StencilInfo.h"

// Hide implementation and avoid collision with GBlk
// TODO: rename cubism_impl::GBlk
namespace cubism_impl {

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
template <class Par_>
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
    static const size_t stridex = Par::bx + 2 * Par::halo;
    static const size_t stridey = Par::by + 2 * Par::halo;
    static const size_t stridez = Par::bz + 2 * Par::halo;
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
        assert(x < stridex * stridey * stridez);

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
class Cubism : public DistrMesh<KF> {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  Cubism(MPI_Comm comm, KF& kf, Vars& par);
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
  using P::par;
  using P::bs_;
  using P::es_; // [fabianw@mavt.ethz.ch; 2019-11-11] not limited anymore (nor power of two)
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

  // hl: number of halo cells from each side
  // cs: number of fields for communication
  static StencilInfo GetStencil(size_t hl, size_t cs) {
    const int a = -int(hl);
    const int b = int(hl) + 1;

    StencilInfo stencil(a,a,a,b,b,b, true, 1, 0);
    stencil.selcomponents.clear();
    for (size_t i = 0; i < cs; ++i) {
      stencil.selcomponents.push_back(i);
    }
    return stencil;
  }
  // Convert Cubism BlockInfo to MyBlockInfo
  static std::vector<MyBlockInfo> GetBlocks(
      const std::vector<BlockInfo>&, MIdx bs, size_t hl);

  std::vector<MIdx> GetBlocks() override;
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

  StreamHdfScal(B& b): b(b) { assert(1 == b.n_comp); }

  // write
  void operate(const int ix, const int iy, const int iz, Scal* out) const
  {
    out[0] = b(ix, iy, iz);
  }

  // read
  void operate(const Scal* out, const int ix, const int iy, const int iz)
  {
    b(ix, iy, iz) = out[0];
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

  StreamHdfVect(B& b): b(b) { assert(3 == b.n_comp); }

  // write
  void operate(const int ix, const int iy, const int iz, Scal* out) const
  {
    const Scal *v = &b(ix, iy, iz);
    out[0] = v[0];
    out[1] = v[1];
    out[2] = v[2];
  }

  // read
  void operate(const Scal* out, const int ix, const int iy, const int iz)
  {
    Scal *v = &b(ix, iy, iz);
    v[0] = out[0];
    v[1] = out[1];
    v[2] = out[2];
  }

  static const char * getAttributeName() { return "Vector"; }
};

// Class with field 'stencil' needed for SynchronizerMPI::sync(Processing)
struct FakeProc {
  StencilInfo stencil;
  explicit FakeProc(StencilInfo si)
    : stencil(si)
  {}
};

template <class Par, class KF>
std::vector<MyBlockInfo> Cubism<Par, KF>::GetBlocks(
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
    b.ptrBlock = c.ptrBlock;
    b.hl = hl;
    bb.push_back(b);
  }
  return bb;
}


template <class Par, class KF>
Cubism<Par, KF>::Cubism(MPI_Comm comm, KF& kf, Vars& par)
  : DistrMesh<KF>(comm, kf, par)
  , g_(p_[0], p_[1], p_[2], b_[0], b_[1], b_[2], ext_, comm)
{
  assert(bs_[0] == FieldView::bx && bs_[1] == FieldView::by &&
      (bs_[2] == FieldView::bz || (bs_[2] == 1 && FieldView::bz == 2)));

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);  // XXX: overwrite isroot_

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Get rid of BlockInfo type
  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // [c]ubism block info
  std::vector<MyBlockInfo> ee = GetBlocks(cc, bs_, hl_);

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
auto Cubism<Par, KF>::GetBlocks() -> std::vector<MIdx> {
  MPI_Barrier(comm_);

  // Get all blocks
  std::vector<BlockInfo> cc = g_.getBlocksInfo(); // all blocks
  auto& m = mk.at(MIdx(cc[0].index))->GetMesh();

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] the comm_id is kernel/comm specific
// and should come from outside
  const std::string comm_id = "test_comm";

  std::vector<BlockInfo> aa;
  // Perform communication if necessary
  if (m.GetComm().size() > 0) {
    // 0. Construct vector of fields to communicate for all blocks on this rank
// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] the fields to be communicated by
// comm_id.  The number and order of fields in this vector should not be changed
// during simulation.  If this vector changes, a new comm_id must be created for
// it.  Ideally, they are created in a std::map outside of GetBlocks()
    std::vector<std::vector<FieldView>> fviews;
    if (g_.isRegistered(comm_id)) {
      fviews = g_.getSynchronizerMPI(comm_id).getFields();
    } else {
      const size_t n_fields = m.GetComm().size();
      fviews.resize(n_fields);
      for (const auto &b : cc) {
        auto &bm = mk.at(MIdx(b.index))->GetMesh();
        auto &bf = bm.GetComm();
        assert(n_fields == bf.size());
        for (size_t fi = 0; fi < n_fields; ++fi) {
          auto &o = bf[fi];
          fviews[fi].push_back(
            FieldView(o->GetBasePtr(), b.index, o->GetSize()));
        }
      }
    }

    size_t cs = 0;
    for (auto& o : m.GetComm()) {
      cs += o->GetSize();
    }

    // 1. Exchange halos in buffer mesh.
    // max(cs, 1) to prevent forbidden call with zero components
    FakeProc fp(GetStencil(hl_, std::max<size_t>(cs, 1)));
    assert(cs == fp.selcomponents.size());

    // Perform communication and register new synchronizer if not already
    // present under label 'comm_id'.  After this point, registered synchronizer
    // objects can be queried with g_.getSynchronizerMPI(comm_id).
    Synch& s = g_.sync(comm_id, fp, fviews);

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

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
template <class Par, class KF>
void Cubism<Par, KF>::ReadBuffer(const std::vector<MIdx>& bb) {}

// FIXME: [fabianw@mavt.ethz.ch; 2019-11-12] Not needed
template <class Par, class KF>
void Cubism<Par, KF>::WriteBuffer(const std::vector<MIdx>& bb) {}

template <class Par, class KF>
void Cubism<Par, KF>::Bcast(const std::vector<MIdx>& bb) {
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
void Cubism<Par, KF>::Scatter(const std::vector<MIdx>& bb) {
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
void Cubism<Par, KF>::Reduce(const std::vector<MIdx>& bb) {
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
void Cubism<Par, KF>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = par.String["dumpformat"];
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
            FieldView(o->GetBasePtr(), b.index, o->GetSize()));
        }
      }

      // Write dump
      auto &bf = m.GetDump();
      assert(bf.size() == fviews.size());
      for (size_t fi = 0; fi < fviews.size(); ++fi) {
        auto fn = GetDumpName(bf[fi].second, "", frame_);
        auto &blocks = fviews[fi];
        if (1 == blocks[0].n_comp) {
          StreamHdfScal<FieldView>::NAME = bf[fi].second;
          DumpHDF5_MPI<StreamHdfScal<FieldView>>(blocks, g_, frame_, frame_, fn,
            ".", Vect(0), m.GetCellSize(), true);
        } else if (3 == blocks[0].n_comp) {
          StreamHdfVect<FieldView>::NAME = bf[fi].second;
          DumpHDF5_MPI<StreamHdfVect<FieldView>>(blocks, g_, frame_, frame_, fn,
            ".", Vect(0), m.GetCellSize(), true);
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
auto Cubism<Par, KF>::GetGlobalBlock() const -> typename M::BlockCells {
  return typename M::BlockCells(p_ * b_ * bs_);
}

template <class Par, class KF>
auto Cubism<Par, KF>::GetGlobalIndex() const -> typename M::IndexCells {
  // TODO: revise, relies on lightweight mesh,
  //       may cause allocation of data fields otherwise
  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx s = GetGlobalBlock().GetSize();
  auto m = InitUniformMesh<M>(dom, MIdx(0), s, 0, true, true, s, 0);
  return m.GetIndexCells();
}

template <class Par, class KF>
auto Cubism<Par, KF>::GetGlobalField(size_t e) -> FieldCell<Scal> {
  using BC = typename M::BlockCells;
  auto gbc = GetGlobalIndex();
  // collective, does actual communication
  auto bb = GetBlocks();
  std::vector<Scal> v(bs_.prod()); // tmp
  MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);
  BC bc(bs_); // cells of one block
  GBlock<size_t, dim> bq(p_ * b_);  // indices of block
  GIndex<size_t, dim> ndq(p_ * b_);  // flat
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
      auto bf = m.GetComm(); // returns vector of shared_ptr
      for (auto& co : m.GetDump()) {
        bf.push_back(co.first);
      }
      assert(e < bf.size());
      const typename M::CoFcs *fcptr = dynamic_cast<typename M::CoFcs*>(bf[e].get());
      assert(fcptr != nullptr);
      const FieldCell<Scal> &fc = *(fcptr->f);
      // copy from inner cells to global field
      for (auto w : bc) {
        gfc[gbc.GetIdx(wb + w)] = fc[mbc.GetIdx(wb + w)];
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
      auto bf = m.GetComm(); // returns vector of shared_ptr
      for (auto& co : m.GetDump()) {
        bf.push_back(co.first);
      }
      assert(e < bf.size());
      const typename M::CoFcs *fcptr = dynamic_cast<typename M::CoFcs*>(bf[e].get());
      assert(fcptr != nullptr);
      const FieldCell<Scal> &fc = *(fcptr->f);
      // copy from inner cells to v
      size_t i = 0;
      // copy from inner cells to global field
      for (auto w : bc) {
        v[i++] = fc[mbc.GetIdx(wb + w)];
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
} // namespace cubism_impl

template <class Scal_, size_t bx_, size_t by_, size_t bz_, size_t es_>
using GPar = cubism_impl::GPar<Scal_, bx_, by_, bz_, es_>;

template <class Par, class KF>
using Cubism = cubism_impl::Cubism<Par, KF>;
