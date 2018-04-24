#pragma once

#include <vector>
#include <limits>
#include <map>
#include <mpi.h>

#include "ILocal.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/output.hpp"
#include "hydro/output_paraview.hpp"
#include "Vars.h"
#include "Distr.h"
#include "hydro/hypre.h"

template <class KF>
class Local : public DistrMesh<KF> {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  using Scal = typename M::Scal;

  Local(MPI_Comm comm, KF& kf, Vars& par);
  typename M::BlockCells GetGlobalBlock() const override;
  // Returns data field i from buffer defined on global mesh
  geom::FieldCell<Scal> GetGlobalField(size_t i) const override; 

 private:
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
  using Rect = geom::Rect<Vect>;
  using IdxCell = geom::IdxCell;
  using P = DistrMesh<KF>;

  using P::mk;
  using P::kf_;
  using P::par;
  using P::bs_;
  using P::es_;
  using P::hl_;
  using P::p_;
  using P::b_; 
  using P::stage_;
  using P::isroot_;
  using P::comm_;
  using P::ext_;

  std::vector<geom::FieldCell<Scal>> buf_; // buffer on mesh
  M gm; // global mesh
  std::unique_ptr<output::Session> session_;
  std::vector<MyBlockInfo> bb_;

  void ReadBuffer(M& m);
  void WriteBuffer(M& m);
  static M CreateMesh(MIdx bs, MIdx b, MIdx p, int es, Scal ext_);
  std::vector<MIdx> GetBlocks() override;
  void ReadBuffer(const std::vector<MIdx>& bb) override;
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void Reduce(const std::vector<MIdx>& bb) override;
  void DumpWrite(const std::vector<MIdx>& bb) override;
};

template <class KF>
auto Local<KF>::CreateMesh(MIdx bs, MIdx b, MIdx p, int es, Scal ext) -> M {
  // Init global mesh
  MIdx ms(bs); // block size 
  MIdx mb(b); // number of blocks
  MIdx mp(p); // number of PEs
  MIdx mm = mp * mb * ms; // total size in cells (without halos)

  Scal h = ext / std::max(std::max(mm[0], mm[1]), mm[2]);
  Vect d0(0); // origin coord
  Vect d1 = d0 + Vect(mm) * h;      // end coord
  Rect d(d0, d1);

  MIdx o(0); // origin index
  std::cout 
    << "o=" << o 
    << " dom=" << d0 << "," << d1 
    << " h=" << h
    << std::endl;

  return geom::InitUniformMesh<M>(d, o, mm, 0);
}

template <class KF>
Local<KF>::Local(MPI_Comm comm, KF& kf, Vars& par) 
  : DistrMesh<KF>(comm, kf, par)
  , buf_(es_)
  , gm(CreateMesh(bs_, b_, p_, es_, ext_))
{

  // Resize buffer for mesh
  for (auto& u : buf_) {
    u.Reinit(gm);
  }

  // Fill block info
  MIdx ms(bs_); // block size 
  MIdx mb(b_[0], b_[1], b_[2]); // number of blocks
  MIdx mp(p_[0], p_[1], p_[2]); // number of PEs
  geom::GBlockCells<3> bc(mb * mp);
  using geom::IdxNode;
  Scal h = (gm.GetNode(IdxNode(1)) - gm.GetNode(IdxNode(0)))[0];
  assert(h > 0);
  std::cerr << "h from gm = " << h << std::endl;
  for (MIdx i : bc) {
    MyBlockInfo b;
    IdxNode n = gm.GetBlockNodes().GetIdx(i * ms);
    Vect o = gm.GetNode(n);
    std::cerr << "o=" << o << " n=" << n.GetRaw() <<  " i=" << i << std::endl;
    for (int q = 0; q < 3; ++q) {
      b.index[q] = i[q];
      b.origin[q] = o[q];
      b.bs[q] = bs_[q];
    }
    b.h_gridpoint = h;
    b.ptrBlock = nullptr;
    b.hl = hl_;
    bb_.push_back(b);
  }

  isroot_ = true; // XXX: overwrite isroot_

  bool islead = true;
  for (auto& b : bb_) {
    MIdx d(b.index);
    b.isroot = (d == MIdx(0));
    b.islead = islead;
    islead = false;
  }

  this->MakeKernels(bb_);
}

template <class KF>
auto Local<KF>::GetBlocks() -> std::vector<MIdx> {
  // Put blocks to map by index 
  std::vector<MIdx> bb;
  for (auto e : bb_) {
    MIdx b(e.index);
    bb.push_back(b);
  }

  return bb;
}

template <class KF>
void Local<KF>::ReadBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();

    ReadBuffer(m);
  }
}

template <class KF>
void Local<KF>::WriteBuffer(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    WriteBuffer(m);
  }
}

template <class KF>
void Local<KF>::Reduce(const std::vector<MIdx>& bb) {
  auto& f = *mk.at(bb[0]); // first kernel
  auto& mf = f.GetMesh();
  auto& vf = mf.GetReduce();  // pointers to reduce

  std::vector<Scal> r(vf.size(), 0); // results

  // Check size is the same for all kernels
  for (auto& b : bb) {
    // TODO: collapse next 2 lines
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  // pointers to reduce
    assert(v.size() == r.size());
  }

  // TODO: Check operation is the same for all kernels

  for (size_t i = 0; i < vf.size(); ++i) {
    Scal r; // result
    std::string s = vf[i].second; // operation string

    enum class Op { sum, prod, max, min };
    Op o;
    if (s == "sum") {
      o = Op::sum;
    } else if (s == "prod") {
      o = Op::prod;
    } else if (s == "max") {
      o = Op::max;
    } else if (s == "min") {
      o = Op::min;
    } else {
      std::cerr << "Reduce(): unknown operation '" << s <<  "'" << std::endl;
      assert(false);
    }

    
    // Init result
    switch (o) {
      case Op::sum:  
        r = 0;  
        break;
      case Op::prod: 
        r = 1;  
        break;
      case Op::max:  
        r = -std::numeric_limits<double>::max();  
        break;
      case Op::min:  
        r = std::numeric_limits<double>::max();  
        break;
      default:
        assert(false);
    }

    // Reduce over all blocks on current rank
    for (auto& b : bb) {
      auto& v = mk.at(b)->GetMesh().GetReduce(); 
      Scal a = *v[i].first;
      switch (o) {
        case Op::sum:  
          r += a;  
          break;
        case Op::prod: 
          r *= a;  
          break;
        case Op::max:  
          r = std::max(r, a);
          break;
        case Op::min:  
          r = std::min(r, a);
          break;
        default:
          assert(false);
      }
    }

    // Write results to all blocks on current rank
    for (auto& b : bb) {
      auto& v = mk.at(b)->GetMesh().GetReduce(); 
      *v[i].first = r;
    }
  }

  // Clear reduce requests
  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    m.ClearReduce();
  }
}

template <class KF>
void Local<KF>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = par.String["dumpformat"];
    if (df == "default") {
      df = "vtk";
    }

    if (df == "vtk") {
      // Initialize on first call
      if (!session_) {
        // TODO: check all blocks are same as first
        output::Content c;
        size_t k = m.GetComm().size() - m.GetDump().size();
        for (auto d : m.GetDump()) {
          c.emplace_back(
              new output::EntryFunction<Scal, IdxCell, M>(
                  d.second, gm, [this,k](IdxCell i) { return buf_[k][i]; }));
          ++k;
        }

        session_.reset(new output::SessionParaviewStructured<M>(
              c, "title", "p" /*filename*/, gm));
      }

      // TODO: Check no change in dump list between time steps
      //       (otherwise session_ needs reinitialization)

      std::cerr << "Dump " << frame_ << ": format=" << df << std::endl;
      ++frame_;
      session_->Write(stage_ * 1., "title"); // TODO: t instead of stage_
    } else {
      P::DumpWrite(bb);
    }
  }
}


template <class KF>
void Local<KF>::ReadBuffer(M& m) {
  int e = 0; // buffer field idx

  for (auto u : m.GetComm()) {
    for (auto i : m.AllCells()) {
      auto& bc = m.GetBlockCells();
      auto& gbc = gm.GetBlockCells();
      MIdx gs = gbc.GetDimensions();
      auto d = bc.GetMIdx(i);
      // periodic
      for (int j = 0; j < 3; ++j) {
        d[j] = (d[j] + gs[j]) % gs[j];
      }
      auto gi = gbc.GetIdx(d);
      (*u)[i] = buf_[e][gi];
    }
    ++e;
  }
}

template <class KF>
void Local<KF>::WriteBuffer(M& m) {
  using MIdx = typename M::MIdx;

  // Check buffer has enough space for all fields
  assert(m.GetComm().size() <= buf_.size() && "Too many fields for Comm()");

  int e = 0; // buffer field idx

  for (auto u : m.GetComm()) {
    for (auto i : m.Cells()) {
      auto& bc = m.GetBlockCells();
      auto& gbc = gm.GetBlockCells();
      auto d = bc.GetMIdx(i); 
      auto gi = gbc.GetIdx(d); 
      buf_[e][gi] = (*u)[i];
    }
    ++e;
  }
}

template <class KF>
auto Local<KF>::GetGlobalBlock() const -> typename M::BlockCells {
  return gm.GetBlockCells();
}

template <class KF>
auto Local<KF>::GetGlobalField(size_t i) const -> geom::FieldCell<Scal> {
  return buf_[i];
}

