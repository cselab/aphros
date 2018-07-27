#pragma once

#include <map>
#include <array>
#include <memory>
#include <mpi.h>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <sstream>

#include "geom/mesh.h"
#include "util/suspender.h"
#include "util/metrics.h"
#include "parse/vars.h"
#include "kernel/kernel.h"
#include "linear/hypre.h"
#include "dump/dump.h"
#include "dump/dumper.h"
#include "report.h"
#include "util/sysinfo.h"

// Abstract block processor.
class Distr {
 public:
  virtual void Run() = 0;
  virtual ~Distr() {}
};

// Abstract block processor aware of Mesh.
// KF: kernel factory derived from KernelMeshFactory
template <class KF>
class DistrMesh : public Distr {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  static constexpr size_t dim = M::dim;
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  virtual void Run();
  virtual void Report();
  virtual ~DistrMesh() {}
  virtual typename M::BlockCells GetGlobalBlock() const;
  virtual typename M::IndexCells GetGlobalIndex() const;
  // Returns data field i from buffer defined on global mesh
  virtual FieldCell<Scal> GetGlobalField(size_t i); 

 protected:
  // TODO: remove comm, needed only by Hypre
  MPI_Comm comm_; // XXX: overwritten by Cubism<KF>
  Vars& par;
  KF& kf_; // kernel factory

  int es_; // element size in Scal
  int hl_; // number of halo cells (same in all directions)
  MIdx bs_; // block size
  MIdx p_; // number of ranks
  MIdx b_; // number of blocks
  Scal ext_; // extent (maximum over all directions)

  int stage_ = 0;
  size_t frame_ = 0; // current dump frame

  bool isroot_; // XXX: overwritten by Local<KF> and Cubism<KF>

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
  // Writes dumps.
  virtual void DumpWrite(const std::vector<MIdx>& bb);
  virtual void ClearComm(const std::vector<MIdx>& bb);
  virtual void ClearDump(const std::vector<MIdx>& bb);
  // TODO: make Pending const
  virtual bool Pending(const std::vector<MIdx>& bb);
  // Create a kernel for each block and put into mk
  // Requires initialized isroot_;
  virtual void MakeKernels(const std::vector<MyBlockInfo>&);
  virtual void TimerReport(const std::vector<MIdx>& bb);
  virtual void ClearTimerReport(const std::vector<MIdx>& bb);

 private:
  MultiTimer<std::string> mt_; // timer all
  MultiTimer<std::string> mtp_; // timer partial
  using LS = typename M::LS;
  std::map<typename LS::T, std::unique_ptr<Hypre>> mhp_; // hypre instances
};

template <class KF>
void DistrMesh<KF>::MakeKernels(const std::vector<MyBlockInfo>& ee) {
  for (auto e : ee) {
    MIdx d(e.index);
    mk.emplace(d, std::unique_ptr<K>(kf_.Make(par, e)));
  }
}

template <class KF>
DistrMesh<KF>::DistrMesh(MPI_Comm comm, KF& kf, Vars& par) 
    : comm_(comm), par(par), kf_(kf)
    , es_(par.Int["comm_size"])
    , hl_(par.Int["hl"])
    , bs_{par.Int["bsx"], par.Int["bsy"], par.Int["bsz"]}
    , p_{par.Int["px"], par.Int["py"], par.Int["pz"]}
    , b_{par.Int["bx"], par.Int["by"], par.Int["bz"]}
    , ext_(par.Double["extent"])
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
void DistrMesh<KF>::TimerReport(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh(); // assume same on all blocks
  std::string fn = m.GetTimerReport();
  if (fn.length()) {
    std::ofstream out;
    out.open(fn);
    out << "mem=" << (sysinfo::GetMem() / double(1 << 20)) << " MB" 
        << std::endl;
    ParseReport(mtp_.GetMap(), out);
    mtp_.Reset();
  }
  ClearTimerReport(bb);
}

template <class KF>
void DistrMesh<KF>::ClearTimerReport(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearTimerReport();
  }
}

template <class KF>
void DistrMesh<KF>::DumpWrite(const std::vector<MIdx>& bb) {
  auto& m = mk.at(bb[0])->GetMesh();
  if (m.GetDump().size()) {
    std::string df = par.String["dumpformat"];
    if (df == "plain") {
      size_t k = 0; // offset in buffer
      // Skip comm 
      for (auto& o : m.GetComm()) {
        k += o->GetSize();
      }
      // Write dump
      for (auto& on : m.GetDump()) {
        std::string fn = GetDumpName(on.second, ".dat", frame_);
        auto ndc = GetGlobalIndex();
        auto bc = GetGlobalBlock();
        auto fc = GetGlobalField(k);
        if (isroot_) {
          Dump(fc, ndc, bc, fn);
        }
        k += on.first->GetSize();
        if (on.first->GetSize() != 1) {
          throw std::runtime_error("DumpWrite(): Support only size 1");
        }
      }
      if (isroot_) {
        std::cerr << "Dump " << frame_ << ": format=" << df << std::endl;
      }
      ++frame_;
    } else {
      throw std::runtime_error("Unknown dumpformat=" + df);
    }
  }
}

// TODO: move
template <class KF>
void DistrMesh<KF>::Solve(const std::vector<MIdx>& bb) {
  const size_t dim = 3;
  auto& vf = mk.at(bb[0])->GetMesh().GetSolve();  // systems to solve on bb[0]

  // Check size is the same for all blocks
  for (auto& b : bb) {
    auto& v = mk.at(b)->GetMesh().GetSolve();  // systems to solve
    if (v.size() != vf.size()) {
      std::stringstream s;
      s << "v.size()=" << v.size() << ",b=" << b
          << " != "
          << "vf.size()=" << vf.size() << ",bf=" << bb[0];
      throw std::runtime_error(s.str());
    }
  }

  for (size_t j = 0; j < vf.size(); ++j) {
    auto& sf = vf[j];  // system to solve on bb[0]
    auto k = sf.t; // key

    if (!mhp_.count(k)) { // create new instance of hypre // XXX
      using LB = typename Hypre::Block;
      std::vector<LB> lbb;
      using LI = typename Hypre::MIdx;

      using MIdx = typename K::MIdx;
      std::vector<MIdx> st = sf.st; // stencil

      for (auto& b : bb) {
        LB lb;
        auto& m = mk.at(b)->GetMesh();
        auto& v = m.GetSolve();
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

      std::string srs = par.String["hypre_symm_solver"]; // solver symm
      assert(srs == "pcg+smg" || srs == "smg" || srs == "pcg" || srs == "zero");

      std::string srg = par.String["hypre_gen_solver"]; // solver gen
      assert(srg == "gmres" || srg == "zero");

      std::string sr; // solver 
      int maxiter;
      using T = typename M::LS::T; // system type
      switch (sf.t) {
        case T::gen:
          sr = srg;
          maxiter = par.Int["hypre_gen_maxiter"];
          break;
        case T::symm:
          sr = srs;
          maxiter = par.Int["hypre_symm_maxiter"];
          break;
        default:
          throw std::runtime_error(
              "Solve(): Unknown system type = " + std::to_string((size_t)sf.t));
      }

      Hypre* hp = new Hypre(comm_, lbb, gs, per, 
                            par.Double["hypre_tol"], par.Int["hypre_print"], 
                            sr, maxiter);
      mhp_.emplace(k, std::unique_ptr<Hypre>(hp));
    } else { // update current instance
      mhp_.at(k)->Update();
    }
    mhp_.at(k)->Solve();
    //mhp_.erase(k);
  }

  for (auto& b : bb) {
    auto& k = *mk.at(b); 
    auto& m = k.GetMesh();
    m.ClearSolve();
  }
}

template <class KF>
bool DistrMesh<KF>::Pending(const std::vector<MIdx>& bb) {
  size_t np = 0;
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    auto& m = k.GetMesh();
    if (m.Pending()) {
      ++np;
    }
  }
  // Check either all blocks done or all pending
  assert(np == 0 || np == bb.size());
  return np;
}

template <class KF>
auto DistrMesh<KF>::GetGlobalBlock() const -> typename M::BlockCells {
  throw std::runtime_error("Not implemented");
  return typename M::BlockCells();
}

template <class KF>
auto DistrMesh<KF>::GetGlobalIndex() const -> typename M::IndexCells {
  throw std::runtime_error("Not implemented");
  return typename M::IndexCells();
}

template <class KF>
auto DistrMesh<KF>::GetGlobalField(size_t) -> FieldCell<Scal> {
  throw std::runtime_error("Not implemented");
  return FieldCell<Scal>();
}

template <class KF>
void DistrMesh<KF>::Run() {
  mt_.Push();
  mtp_.Push();
  do {
    auto bb = GetBlocks();
    
    assert(!bb.empty());

    ReadBuffer(bb);

    DumpWrite(bb);
    ClearDump(bb);

    ClearComm(bb);

    Reduce(bb);

    Solve(bb);

    mt_.Pop(mk.at(bb[0])->GetMesh().GetCurName());
    mt_.Push();

    mtp_.Pop(mk.at(bb[0])->GetMesh().GetCurName());
    TimerReport(bb);
    mtp_.Push();

    Run(bb);

    // Write comm and dump
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
  mt_.Pop("last");
  mtp_.Pop("last");

  Report();
}

template <class KF>
void DistrMesh<KF>::Report() {
  if (isroot_) {
    double a = 0.; // total
    for (auto e : mt_.GetMap()) {
      a += e.second;
    }

    if (par.Int["verbose_stages"]) {
      std::cout << std::fixed;
      auto& mp = mt_.GetMap();
      std::cout 
          << "mem=" << (sysinfo::GetMem() / double(1 << 20)) << " MB" 
          << std::endl;
      ParseReport(mp, std::cout);
      for (auto e : mp) {
        auto n = e.first; // name
        if (n == "") {
          n = "other";
        }
        auto t = e.second; // time

        std::cout 
            << n << "\n" 
            << std::setprecision(5) << t << " s = "
            << std::setprecision(3) << 100. * t / a << "%\n";
      }
      std::cout << std::endl;
    }

    MIdx gs;
    {
      MIdx p(par.Int["px"], par.Int["py"], par.Int["pz"]);
      MIdx b(par.Int["bx"], par.Int["by"], par.Int["bz"]);
      MIdx bs(par.Int["bsx"], par.Int["bsy"], par.Int["bsz"]);
      gs = p * b * bs;
    }
    size_t nc = gs.prod(); // total cells
    size_t nt = par.Int["max_step"];
    size_t ni = par.Int["iter"];

    // Returns: hour, minute, second, millisecond
    auto get_hmsm = [](double t) {
      std::array<int, 4> r;
      r[0] = int(t / 3600);
      t -= r[0] * 3600;
      r[1] = int(t / 60);
      t -= r[1] * 60;
      r[2] = int(t);
      r[3] = int((t - int(t)) * 1000) % 1000;
      return r;
    };

    auto h = get_hmsm(a);
    std::cout << std::setprecision(5) << std::scientific
        << "cells = " << nc << "\n"
        << "steps = " << nt << "\n"
        << "iters = " << ni << "\n"
        << "total = " << int(a) << " s"
        << " = " << std::setfill('0')
        << std::setw(2) << h[0] << ":" 
        << std::setw(2) << h[1] << ":" 
        << std::setw(2) << h[2] << "." 
        << std::setw(3) << h[3] << "\n"
        << "time/cell/iter = " << a / (nc * ni) << " s"
        << std::endl;
  }
}
