#pragma once

#include <vector>
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

using Scal = double;

template <class KF>
class Local : public DistrMesh<KF> {
 public:
  using K = typename KF::K;
  using M = typename KF::M;

  Local(MPI_Comm comm, KF& kf, Vars& par);
  typename M::BlockCells GetGlobalBlock() const override;
  // Returns data field i from buffer defined on global mesh
  geom::FieldCell<Scal> GetGlobalField(size_t i) const override; 

 private:
  using MIdx = typename  M::MIdx;
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
  using P:: b_; 
  using P::step_;
  using P::stage_;
  using P::frame_;
  using P::isroot_;
  using P::comm_;

  std::vector<geom::FieldCell<Scal>> buf_; // buffer on mesh
  M gm; // global mesh
  std::unique_ptr<output::Session> session_;
  std::vector<MyBlockInfo> bb_;

  void ReadBuffer(M& m);
  void WriteBuffer(M& m);
  static M CreateMesh(MIdx bs, MIdx b, MIdx p, int es);
  std::vector<MIdx> GetBlocks() override;
  void ReadBuffer(const std::vector<MIdx>& bb) override;
  void WriteBuffer(const std::vector<MIdx>& bb) override;
  void Reduce(const std::vector<MIdx>& bb) override;
  void Dump(int frame, int step) override;
};

template <class KF>
auto Local<KF>::CreateMesh(MIdx bs, MIdx b, MIdx p, int es) -> M {
  // Init global mesh
  MIdx ms(bs); // block size 
  MIdx mb(b); // number of blocks
  MIdx mp(p); // number of PEs
  MIdx mm = mp * mb * ms; // total size in cells (without halos)

  Scal ext = 1.; // TODO: extent from par
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
  , gm(CreateMesh(bs_, b_, p_, es_))
{

  output::Content content = {
    std::make_shared<output::EntryFunction<Scal, IdxCell, M>>(
        "vx", gm, [this](IdxCell i) { return buf_[0][i]; }),
    std::make_shared<output::EntryFunction<Scal, IdxCell, M>>(
        "vy", gm, [this](IdxCell i) { return buf_[1][i]; }),
    std::make_shared<output::EntryFunction<Scal, IdxCell, M>>(
        "vz", gm, [this](IdxCell i) { return buf_[2][i]; }),
    std::make_shared<output::EntryFunction<Scal, IdxCell, M>>(
        "p", gm, [this](IdxCell i) { return buf_[3][i]; })
  };

  session_.reset(new output::SessionParaviewStructured<M>(
        content, "title", "p" /*filename*/, gm));

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
    }
    b.h_gridpoint = h;
    b.ptrBlock = nullptr;
    b.hl = hl_;
    bb_.push_back(b);
  }

  for (auto e : bb_) {
    auto d = e.index;
    // TODO: constructor
    MIdx b(d[0], d[1], d[2]);
    mk.emplace(b, std::unique_ptr<K>(kf_.Make(par, e)));
  }

  isroot_ = true;
}

template <class KF>
auto Local<KF>::GetBlocks() -> std::vector<MIdx> {
  // Put blocks to map by index 
  std::vector<MIdx> bb;
  for (auto e : bb_) {
    auto d = e.index;
    // TODO: constructor
    MIdx b(d[0], d[1], d[2]);
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

  // Reduce over all kernels on current rank
  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  
    for (size_t i = 0; i < r.size(); ++i) {
      r[i] += *v[i];
    }
  }

  // Write results to all kernels on current rank
  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    auto& v = m.GetReduce();  
    for (size_t i = 0; i < r.size(); ++i) {
      *v[i] = r[i];
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
void Local<KF>::Dump(int frame, int step) {
  auto suff = "_" + std::to_string(frame);
  std::cerr << "Output" << std::endl;
  session_->Write(step * 1., "title:0");
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

  m.ClearComm();
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
geom::FieldCell<Scal> Local<KF>::GetGlobalField(size_t i) const {
  return buf_[i];
}

