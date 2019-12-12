#pragma once

#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <array>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "dump/dump.h"
#include "dump/dumper.h"
#include "geom/mesh.h"
#include "linear/hypre.h"
#include "linear/hypresub.h"
#include "parse/vars.h"
#include "report.h"
#include "util/histogram.h"
#include "util/metrics.h"
#include "util/suspender.h"
#include "util/sysinfo.h"

// Abstract block processor aware of Mesh.
// KF: kernel factory derived from KernelMeshFactory
template <class KF>
class DistrMesh {
 public:
  using K = typename KF::K;
  using M = typename KF::M;
  static constexpr size_t dim = M::dim;
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  virtual void Run();
  virtual void Report();
  virtual void ReportOpenmp();
  virtual ~DistrMesh() {
    hist_.Append(this->samp_); // collect samples from derived classes
    const auto it_comp = hist_.GetSamples().find("RunKernels(inner)");
    if (it_comp != hist_.GetSamples().end()) {
      hist_.Insert("ComputeTime", it_comp->second);
    }
    const auto it_wait = hist_.GetSamples().find("waitall_avail_halo");
    if (it_wait != hist_.GetSamples().end()) {
      hist_.Insert("TransferTime", it_wait->second);
    }
  }
  virtual typename M::BlockCells GetGlobalBlock() const;
  virtual typename M::IndexCells GetGlobalIndex() const;
  // Returns data field i from buffer defined on global mesh
  virtual FieldCell<Scal> GetGlobalField(size_t i);

 protected:
  // TODO: remove comm, needed only by Hypre
  MPI_Comm comm_; // XXX: overwritten by Cubism<KF>
  const Vars& var;
  Vars& var_mutable;
  const KF& kf_; // kernel factory
  Sampler samp_; // sampler accessible to derived classes

  int hl_; // number of halo cells (same in all directions)
  MIdx bs_; // block size
  MIdx p_; // number of ranks
  MIdx b_; // number of blocks
  Scal ext_; // extent (maximum over all directions)

  int stage_ = 0;
  size_t frame_ = 0; // current dump frame

  bool isroot_; // XXX: overwritten by Local<KF> and Cubism<KF>

  std::map<MIdx, std::unique_ptr<K>, typename MIdx::LexLess> mk;

  DistrMesh(MPI_Comm comm, KF& kf, Vars& var);
  // Performs communication and returns indices of blocks with updated halos.
  virtual std::vector<MIdx> GetBlocks(bool inner) = 0;
  virtual std::vector<MIdx> GetBlocks() {
    auto bbi = GetBlocks(true);
    auto bbh = GetBlocks(false);
    bbi.insert(bbi.end(), bbh.begin(), bbh.end());
    return bbi;
  }
  // Fill selected halo cells with garbage
  void ApplyNanFaces(const std::vector<MIdx>& bb);
  // Copy from communication buffer to fields
  virtual void ReadBuffer(const std::vector<MIdx>& bb) = 0;
  // Call kernels for current stage
  virtual void RunKernels(const std::vector<MIdx>& bb);
  // Copy from fields to communication buffer
  virtual void WriteBuffer(const std::vector<MIdx>& bb) = 0;
  // Reduce TODO: extend doc
  virtual void Reduce(const std::vector<MIdx>& bb) = 0;
  virtual void Scatter(const std::vector<MIdx>& bb) = 0;
  virtual void Bcast(const std::vector<MIdx>& bb) = 0;
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
  Histogram hist_; // histogram sample collector
  MultiTimer<std::string> mt_; // timer all
  MultiTimer<std::string> mtp_; // timer partial
  using LS = typename M::LS;
  std::map<typename LS::T, std::unique_ptr<HypreSub>> mhp_; // hypre instances
};

template <class KF>
void DistrMesh<KF>::MakeKernels(const std::vector<MyBlockInfo>& ee) {
  for (auto e : ee) {
    MIdx d(e.index);
    mk.emplace(d, std::unique_ptr<K>(kf_.Make(var_mutable, e)));
  }
}

template <class KF>
DistrMesh<KF>::DistrMesh(MPI_Comm comm, KF& kf, Vars& var0)
    : comm_(comm)
    , var(var0)
    , var_mutable(var0)
    , kf_(kf)
    , samp_(var.Int["histogram"])
    , hl_(var.Int["hl"])
    , bs_{var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]}
    , p_{var.Int["px"], var.Int["py"], var.Int["pz"]}
    , b_{var.Int["bx"], var.Int["by"], var.Int["bz"]}
    , ext_(var.Double["extent"])
    , hist_(comm, "distrmesh", var.Int["histogram"]) {}

template <class KF>
void DistrMesh<KF>::RunKernels(const std::vector<MIdx>& bb) {
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < bb.size(); ++i) {
    mk.at(bb[i])->Run();
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
    std::string df = var.String["dumpformat"];
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

      LI per{false, false, false};
      per[0] = var.Int["hypre_periodic_x"];
      per[1] = var.Int["hypre_periodic_y"];
      per[2] = var.Int["hypre_periodic_z"];

      LI gs(bs_ * b_ * p_);

      mhp_.emplace(k, new HypreSub(comm_, lbb, gs, per));
    } else { // update current instance
      mhp_.at(k)->Update();
    }

    auto& s = mhp_.at(k);

    std::string sr; // solver
    int maxiter;
    Scal tol;

    {
      std::string srs = var.String["hypre_symm_solver"]; // solver symm
      assert(srs == "pcg+smg" || srs == "smg" || srs == "pcg" || srs == "zero");

      std::string srg = var.String["hypre_gen_solver"]; // solver gen
      assert(srg == "gmres" || srg == "zero");

      using T = typename M::LS::T; // system type
      switch (sf.t) {
        case T::gen:
          sr = srg;
          maxiter = var.Int["hypre_gen_maxiter"];
          tol = var.Double["hypre_gen_tol"];
          break;
        case T::symm:
          sr = srs;
          maxiter = var.Int["hypre_symm_maxiter"];
          tol = var.Double["hypre_symm_tol"];
          break;
        default:
          throw std::runtime_error(
              "Solve(): Unknown system type = " + std::to_string((size_t)sf.t));
      }

    }

    if (sf.prefix != "") {
      sr = var.String["hypre_" + sf.prefix + "_solver"];
      maxiter = var.Int["hypre_" + sf.prefix + "_maxiter"];
      tol = var.Double["hypre_" + sf.prefix + "_tol"];
    }

    hist_.SeedSample();
    s->Solve(tol, var.Int["hypre_print"], sr, maxiter);
    hist_.CollectSample("Hypre::Solve");

    for (auto& b : bb) {
      auto& m = mk.at(b)->GetMesh();
      m.SetResidual(s->GetResidual());
      m.SetIter(s->GetIter());
    }
  }

  for (auto& b : bb) {
    auto& m = mk.at(b)->GetMesh();
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
void DistrMesh<KF>::ApplyNanFaces(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = mk.at(b)->GetMesh();
    for (auto& o : m.GetComm()) {
      if (auto od = dynamic_cast<typename M::CoFcs*>(o.get())) {
        m.ApplyNanFaces(*od->f);
      } else if (auto od = dynamic_cast<typename M::CoFcv*>(o.get())) {
        m.ApplyNanFaces(*od->f);
      } else {
        throw std::runtime_error("Distr::Run(): unknown field type for nan");
      }
    }
  }
}

template <class KF>
void DistrMesh<KF>::Run() {
  if (var.Int["verbose_openmp"]) {
    ReportOpenmp();
  }
  mt_.Push();
  mtp_.Push();
  do {
    hist_.SeedSample();
    std::vector<MIdx> bb;
    if (mk.begin()->second->GetMesh().GetDump().size() > 0) {
      hist_.SeedSample();
      bb = GetBlocks(); // all blocks, sync communication
      hist_.CollectSample("GetBlocks(sync)");
      ReadBuffer(bb);
      ApplyNanFaces(bb);
      DumpWrite(bb);
      ClearDump(bb);
      ClearComm(bb);
      hist_.SeedSample();
      RunKernels(bb);
      hist_.CollectSample("RunKernels(sync)");
    } else {
      hist_.SeedSample();
      auto bbi = GetBlocks(true); // inner blocks, async communication
      hist_.CollectSample("GetBlocks(inner)");
      ReadBuffer(bbi);
      ApplyNanFaces(bbi);
      ClearComm(bbi);
      hist_.SeedSample();
      RunKernels(bbi);
      hist_.CollectSample("RunKernels(inner)");

      hist_.SeedSample();
      auto bbh = GetBlocks(false); // halo blocks, wait for communication
      hist_.CollectSample("GetBlocks(halo)");
      ReadBuffer(bbh);
      ApplyNanFaces(bbh);
      ClearComm(bbh);
      hist_.SeedSample();
      RunKernels(bbh);
      hist_.CollectSample("RunKernels(halo)");

      bb = bbi;
      bb.insert(bb.end(), bbh.begin(), bbh.end());
    }

    WriteBuffer(bb);

    stage_ += 1;

    // Print current stage name
    if (isroot_ && var.Int["verbose"]) {
      auto& m = mk.begin()->second->GetMesh();
      std::cerr << "*** STAGE"
          << " #" << stage_
          << " depth=" << m.GetDepth()
          << " " << m.GetCurName()
          << " ***" << std::endl;
    }

    // Break if no pending stages
    if (!Pending(bb)) {
      hist_.CollectSample("Run");
      hist_.PopLast("Run"); // last sample is invalid
      break;
    }

    hist_.SeedSample();
    Reduce(bb);
    hist_.CollectSample("Reduce");

    hist_.SeedSample();
    Scatter(bb);
    hist_.CollectSample("Scatter");

    hist_.SeedSample();
    Bcast(bb);
    hist_.CollectSample("Bcast");

    Solve(bb);

    hist_.CollectSample("Run");

    mt_.Pop(mk.at(bb[0])->GetMesh().GetCurName());
    mt_.Push();

    mtp_.Pop(mk.at(bb[0])->GetMesh().GetCurName());
    TimerReport(bb);
    mtp_.Push();
  } while (true);
  mt_.Pop("last");
  mtp_.Pop("last");

  std::vector<MIdx> bb = GetBlocks();
  for (const auto& b : bb) {
    const auto& samp = mk.at(b)->GetMesh().GetSampler();
    hist_.Append(samp); // collect samples from meshes
  }

  if (var.Int["verbose_time"]) {
    Report();
  }
}

template <class KF>
void DistrMesh<KF>::Report() {
  if (isroot_) {
    double a = 0.; // total
    for (auto e : mt_.GetMap()) {
      a += e.second;
    }

    if (var.Int["verbose_stages"]) {
      std::cout << std::fixed;
      auto& mp = mt_.GetMap();
      std::cout
          << "mem=" << (sysinfo::GetMem() / double(1 << 20)) << " MB"
          << std::endl;
      ParseReport(mp, std::cout);
    }

    MIdx gs;
    {
      MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
      MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
      MIdx bs(var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]);
      gs = p * b * bs;
    }
    size_t nc = gs.prod(); // total cells
    size_t nt = var.Int["max_step"];
    size_t ni = var.Int["iter"];

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

template <class KF>
void DistrMesh<KF>::ReportOpenmp() {
  #ifdef _OPENMP
  if (isroot_) {
    #pragma omp parallel
    {
      #pragma omp single
      std::cout << "OpenMP threads" << std::endl;
      #pragma omp for ordered
      for (int i = 0; i < omp_get_num_threads(); ++i) {
        #pragma omp ordered
        {
          std::cout << "thread="
                    << std::setw(2) << omp_get_thread_num();
          std::cout << std::setw(8) << " cpu="
                    << std::setw(2) << sched_getcpu();
          std::cout << std::endl;
        }
      }
    }
  }
  #endif
}
