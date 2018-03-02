#pragma once

#include <map>
#include <array>
#include <memory>
#include <mpi.h>

#include "Vars.h"
#include "Kernel.h"
#include "hydro/hypre.h"
#include "hydro/suspender.h"


class Distr {
 public:
  virtual bool IsDone() const = 0;
  virtual void Step() = 0;
  virtual ~Distr() {}
};

template <class KF>
class DistrMesh : public Distr {
 public:
  virtual bool IsDone() const;
  virtual void Step();
  virtual ~DistrMesh() {}

 protected:
  using K = typename KF::K;
  using M = typename KF::M;
  static constexpr size_t dim = M::dim;
  using MIdx = typename M::MIdx;

  // TODO: remove, needed only by Hypre
  // XXX: overwritten by Cubism<KF>
  MPI_Comm comm_;
  Vars& par;
  KF& kf_; // kernel factory

  int bs_; // block size
  int es_; // element size in Scal
  int hl_; // number of halo cells (same in all directions)
  MIdx p_; // number of ranks
  MIdx b_; // number of blocks

  int step_ = 0;
  int stage_ = 0;
  int frame_ = 0;

  // XXX: overwritten by Local<KF>
  bool isroot_;

  std::map<MIdx, std::unique_ptr<K>, typename MIdx::LexLess> mk;

  DistrMesh(MPI_Comm comm, KF& kf, int bs, int es, int hl, Vars& par);
  virtual std::vector<MIdx> GetBlocks() = 0;
  // Copy data from buffer halos to fields collected by Comm()
  virtual void ReadBuffer(const std::vector<MIdx>& bb) = 0;
  // Call kernels for current stage
  virtual void Run(const std::vector<MIdx>& bb);
  // Copy data to buffer mesh from fields collected by Comm()
  virtual void WriteBuffer(const std::vector<MIdx>& bb) = 0;
  // Reduce TODO: extend
  virtual void Reduce(const std::vector<MIdx>& bb) = 0;
  // Solve TODO: extend
  virtual void Solve(const std::vector<MIdx>& bb);
  virtual bool Pending(const std::vector<MIdx>& bb);
  virtual void Dump(int frame, int step) = 0;
};

template <class KF>
DistrMesh<KF>::DistrMesh(MPI_Comm comm, KF& kf, int bs, int es, int hl, Vars& par) 
  : comm_(comm), par(par), kf_(kf), bs_(bs), es_(es), hl_(hl)
  , p_{par.Int["px"], par.Int["py"], par.Int["pz"]}
  , b_{par.Int["bx"], par.Int["by"], par.Int["bz"]}
{}

template <class KF>
bool DistrMesh<KF>::IsDone() const { 
  return step_ > par.Int["max_step"]; 
}

template <class KF>
void DistrMesh<KF>::Run(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    k.Run();
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
      using B = MyBlock;
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
    if (par.Int["hypre_periodic"]) {
      for (size_t i = 0; i < dim; ++i) {
        per[i] = true;
      }
    }

    LI gs(dim);
    for (size_t i = 0; i < dim; ++i) {
      gs[i] = bs_ * b_[i] * p_[i];
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
void DistrMesh<KF>::Step() {
  if (isroot_) {
    std::cerr << "***** STEP " << step_ << " ******" << std::endl;
  }
  stage_ = 0;
  do {
    auto bb = GetBlocks();
    
    assert(!bb.empty());
  
    ReadBuffer(bb);
    
    Run(bb);

    WriteBuffer(bb);

    Reduce(bb);

    Solve(bb);

    stage_ += 1;

    // Print current stage name
    if (isroot_) {
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

  if (step_ % (par.Int["max_step"] / par.Int["num_frames"])  == 0) {
    Dump(frame_, step_);
    ++frame_;
  }
  ++step_;
}

