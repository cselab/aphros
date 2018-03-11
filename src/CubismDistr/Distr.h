#pragma once

#include <map>
#include <array>
#include <memory>
#include <mpi.h>

#include "Vars.h"
#include "Kernel.h"
#include "hydro/hypre.h"
#include "hydro/suspender.h"
#include "hydro/mesh.hpp"


class Distr {
 public:
  virtual void Run() = 0;
  virtual ~Distr() {}
};

template <class KF>
class DistrMesh : public Distr {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  static constexpr size_t dim = M::dim;
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;

  virtual void Run();
  virtual ~DistrMesh() {}
  virtual typename M::BlockCells GetGlobalBlock() const;
  // Returns data field i from buffer defined on global mesh
  virtual geom::FieldCell<Scal> GetGlobalField(size_t i) const; 

 protected:

  // TODO: remove, needed only by Hypre
  // XXX: overwritten by Cubism<KF>
  MPI_Comm comm_;
  Vars& par;
  KF& kf_; // kernel factory

  MIdx bs_; // block size
  int es_; // element size in Scal
  int hl_; // number of halo cells (same in all directions)
  MIdx p_; // number of ranks
  MIdx b_; // number of blocks

  int stage_ = 0;
  int frame_ = 0;

  // XXX: overwritten by Local<KF>
  bool isroot_;

  std::map<MIdx, std::unique_ptr<K>, typename MIdx::LexLess> mk;

  DistrMesh(MPI_Comm comm, KF& kf, Vars& par);
  virtual std::vector<MIdx> GetBlocks() = 0;
  // Copy data from buffer halos to fields collected by Comm()
  virtual void ReadBuffer(const std::vector<MIdx>& bb) = 0;
  // Call kernels for current stage
  virtual void Run(const std::vector<MIdx>& bb);
  // Copy data to buffer mesh from fields collected by Comm()
  virtual void WriteBuffer(const std::vector<MIdx>& bb) = 0;
  // Reduce TODO: extend doc
  virtual void Reduce(const std::vector<MIdx>& bb) = 0;
  virtual void Solve(const std::vector<MIdx>& bb);
  virtual void DumpComm(const std::vector<MIdx>& bb);
  virtual void DumpWrite(const std::vector<MIdx>& bb) = 0;
  virtual void ClearComm(const std::vector<MIdx>& bb);
  virtual void ClearDump(const std::vector<MIdx>& bb);
  // TODO: make Pending const
  virtual bool Pending(const std::vector<MIdx>& bb);
};

template <class KF>
DistrMesh<KF>::DistrMesh(MPI_Comm comm, KF& kf, Vars& par) 
  : comm_(comm), par(par), kf_(kf)
  , es_(par.Int["comm_size"])
  , hl_(par.Int["hl"])
  , bs_{par.Int["bsx"], par.Int["bsy"], par.Int["bsz"]}
  , p_{par.Int["px"], par.Int["py"], par.Int["pz"]}
  , b_{par.Int["bx"], par.Int["by"], par.Int["bz"]}
{}

template <class KF>
void DistrMesh<KF>::Run(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    k.Run();
    }
}

template <class KF>
void DistrMesh<KF>::ClearComm(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearComm();
  }
}

template <class KF>
void DistrMesh<KF>::ClearDump(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearDump();
  }
}

template <class KF>
void DistrMesh<KF>::DumpComm(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = mk.at(b)->GetMesh();

    for (auto d : m.GetDump()) {
      m.Comm(d.first);
    }
  }
}


template <class KF>
void DistrMesh<KF>::Solve(const std::vector<MIdx>& bb) {
  const size_t dim = 3;
  auto& f = *mk.at(bb[0]); // first kernel
  auto& mf = f.GetMesh();
  auto& vf = mf.GetSolve();  // LS to solve

  // Check size is the same for all kernels
  for (auto& b : bb) {
    auto& k = *mk.at(b); // kernel
    auto& m = k.GetMesh();
    auto& v = m.GetSolve();  // pointers to reduce
    assert(v.size() == vf.size());
  }

  for (size_t j = 0; j < vf.size(); ++j) {
    using LB = typename Hypre::Block;
    std::vector<LB> lbb;
    using LI = typename Hypre::MIdx;

    using MIdx = typename K::MIdx;
    std::vector<MIdx> st = vf[j].st; // stencil

    for (auto& b : bb) {
      LB lb;
      auto& k = *mk.at(b); // kernel
      auto& m = k.GetMesh();
      auto& v = m.GetSolve();  // pointers to reduce
      auto& bc = m.GetInBlockCells();
      auto& s = v[j];  
      lb.l = bc.GetBegin();
      lb.u = bc.GetEnd() - MIdx(1);
      for (MIdx& e : st) {
        lb.st.emplace_back(e);
      }
      lb.a = s.a;
      lb.r = s.b;
      lb.x = s.x;
      lbb.push_back(lb);
    }

    std::vector<bool> per(dim, false);
    per[0] = par.Int["hypre_periodic_x"];
    per[1] = par.Int["hypre_periodic_y"];
    per[2] = par.Int["hypre_periodic_z"];

    LI gs(dim);
    for (size_t i = 0; i < dim; ++i) {
      gs[i] = bs_[i] * b_[i] * p_[i];
    }

    Hypre hp(comm_, lbb, gs, per, 
             par.Double["hypre_tol"], par.Int["hypre_print"]);
    hp.Solve();
  }

  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    m.ClearSolve();
  }
}

template <class KF>
bool DistrMesh<KF>::Pending(const std::vector<MIdx>& bb) {
  int np = 0;
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    auto& m = k.GetMesh();
    if (m.Pending()) {
      ++np;
    }
  }
  // Check either all done or all pending
  assert(np == 0 || np == bb.size());
  return np;
}

template <class KF>
auto DistrMesh<KF>::GetGlobalBlock() const -> typename M::BlockCells {
  assert(false && "Not implemented");
  return typename M::BlockCells();
}

template <class KF>
geom::FieldCell<Scal> DistrMesh<KF>::GetGlobalField(size_t i) const {
  assert(false && "Not implemented");
  return geom::FieldCell<Scal>();
}


template <class KF>
void DistrMesh<KF>::Run() {
  stage_ = 0;
  do {
    auto bb = GetBlocks();
    
    assert(!bb.empty());

    ReadBuffer(bb);

    DumpWrite(bb);
    ClearDump(bb);

    ClearComm(bb);

    Reduce(bb);

    Solve(bb);

    Run(bb);

    DumpComm(bb);

    WriteBuffer(bb);

    stage_ += 1;

    // Print current stage name
    if (isroot_ && par.Int["verbose"]) {
      auto& m = mk.at(bb[0])->GetMesh();
      std::cerr << "*** STAGE"
          << " #" << stage_ 
          << " depth=" << m.GetDepth() 
          << " " << m.GetCurName() 
          << " ***" << std::endl;
    }

    // Break if no pending stages
    if (!Pending(bb)) {
      break;
    }
  } while (true);

  /*
  if (par.Int["output"] && 
      step_ % (par.Int["max_step"] / par.Int["num_frames"])  == 0) {
    Dump(frame_, step_);
    ++frame_;
  }
  */
}

