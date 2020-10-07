// Created by Petr Karnakov on 13.12.2019
// Copyright 2019 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <stdexcept>

#include "distr.h"

template <class M>
void DistrMesh<M>::MakeKernels(const std::vector<MyBlockInfo>& ee) {
  for (auto e : ee) {
    const MIdx d(e.index);
    std::unique_ptr<KernelMesh<M>> kernel(kf_.Make(var_mutable, e));
    auto& m = kernel->GetMesh();
    m.flags.comm = comm_;
    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    m.flags.is_periodic[2] = var.Int["hypre_periodic_z"];
    mk.emplace(d, kernel.release());
  }
}

template <class M>
DistrMesh<M>::DistrMesh(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var0)
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

template <class M>
DistrMesh<M>::~DistrMesh() {
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

template <class M>
void DistrMesh<M>::RunKernels(const std::vector<MIdx>& bb) {
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < bb.size(); ++i) {
    try {
      mk.at(bb[i])->Run();
    } catch (const std::exception& e) {
      std::cerr << "\nabort on block=" << i << " after throwing exception\n"
                << e.what() << '\n';
      throw e;
    }
  }
}

template <class M>
void DistrMesh<M>::ClearComm(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearComm();
  }
}

template <class M>
void DistrMesh<M>::ClearDump(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearDump();
  }
}

template <class M>
void DistrMesh<M>::TimerReport(const std::vector<MIdx>& bb) {
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

template <class M>
void DistrMesh<M>::ClearTimerReport(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearTimerReport();
  }
}

template <class M>
void DistrMesh<M>::ReduceToLead(const std::vector<MIdx>& bb) {
  using ReductionScal = typename M::OpS;
  using ReductionScalInt = typename M::OpSI;
  using ReductionConcat = typename M::OpCat;
  auto& first = *mk.at(bb[0]); // first kernel
  auto& mfirst = first.GetMesh();
  auto& vfirst = mfirst.GetReduceToLead(); // reduce requests

  // Check size is the same for all kernels
  for (auto& b : bb) {
    const auto& v = mk.at(b)->GetMesh().GetReduceToLead(); // reduce requests
    fassert_equal(v.size(), vfirst.size());
  }

  auto append = [&bb, this](auto* derived, auto& buf, size_t i) {
    for (auto& b : bb) {
      const auto& v = mk.at(b)->GetMesh().GetReduceToLead();
      dynamic_cast<decltype(derived)>(v[i].get())->Append(buf);
    }
  };
  auto set = [&bb, this](auto* derived, auto& buf, size_t i) {
    for (auto& b : bb) {
      const auto& m = mk.at(b)->GetMesh();
      if (m.IsLead()) {
        const auto& v = m.GetReduceToLead();
        dynamic_cast<decltype(derived)>(v[i].get())->Set(buf);
      }
    }
  };
  auto apply = [append, set](auto* derived, size_t i) {
    if (derived) {
      auto buf = derived->Neutral();
      append(derived, buf, i);
      set(derived, buf, i);
      return true;
    }
    return false;
  };

  for (size_t i = 0; i < vfirst.size(); ++i) {
    auto* request = vfirst[i].get();
    if (apply(dynamic_cast<ReductionScal*>(request), i) ||
        apply(dynamic_cast<ReductionScalInt*>(request), i) ||
        apply(dynamic_cast<ReductionConcat*>(request), i)) {
    } else {
      throw std::runtime_error(
          FILELINE + "ReduceToLead: Unknown M::Op implementation");
    }
  }

  // Clear reduce requests
  for (auto& b : bb) {
    mk.at(b)->GetMesh().ClearReduceToLead();
  }
}

template <class M>
void DistrMesh<M>::DumpWrite(const std::vector<MIdx>& bb) {
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

template <class M>
bool DistrMesh<M>::Pending(const std::vector<MIdx>& bb) {
  size_t np = 0;
  for (auto& b : bb) {
    auto& k = *mk.at(b);
    auto& m = k.GetMesh();
    if (m.Pending()) {
      ++np;
    }
  }
  // Check either all blocks done or all pending
  fassert(np == 0 || np == bb.size());
  return np;
}

template <class M>
auto DistrMesh<M>::GetGlobalBlock() const -> typename M::BlockCells {
  throw std::runtime_error("Not implemented");
  return typename M::BlockCells();
}

template <class M>
auto DistrMesh<M>::GetGlobalIndex() const -> typename M::IndexCells {
  throw std::runtime_error("Not implemented");
  return typename M::IndexCells();
}

template <class M>
auto DistrMesh<M>::GetGlobalField(size_t) -> FieldCell<Scal> {
  throw std::runtime_error("Not implemented");
  return FieldCell<Scal>();
}

template <class M>
void DistrMesh<M>::ApplyNanFaces(const std::vector<MIdx>& bb) {
  for (auto& b : bb) {
    auto& m = mk.at(b)->GetMesh();
    for (auto& o : m.GetComm()) {
      if (auto* os = dynamic_cast<typename M::CoFcs*>(o.get())) {
        m.ApplyNanFaces(*os->field);
      } else if (auto* ov = dynamic_cast<typename M::CoFcv*>(o.get())) {
        m.ApplyNanFaces(*ov->field);
      } else {
        throw std::runtime_error("Distr::Run(): unknown field type for nan");
      }
    }
  }
}

template <class M>
void DistrMesh<M>::Run() {
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
      const auto& mf = mk.begin()->second->GetMesh();
      std::cerr << "*** STAGE"
                << " #" << stage_ << " depth=" << mf.GetDepth() << " "
                << mf.GetCurName() << " ***" << std::endl;
      // Check stage name is the same for all kernels
      for (auto& b : bb) {
        auto& m = mk.at(b)->GetMesh();
        if (m.GetCurName() != mf.GetCurName()) {
          std::stringstream s;
          s << "Blocks " << b << " and " << mk.begin()->first
            << " diverged to different stages '" << m.GetCurName() << "' and '"
            << mf.GetCurName() << "'";
          throw std::runtime_error(s.str());
        }
      }
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
    ReduceToLead(bb);
    hist_.CollectSample("ReduceToLead");

    hist_.SeedSample();
    Scatter(bb);
    hist_.CollectSample("Scatter");

    hist_.SeedSample();
    Bcast(bb);
    hist_.CollectSample("Bcast");

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

template <class M>
void DistrMesh<M>::Report() {
  if (isroot_) {
    double a = 0.; // total
    for (auto e : mt_.GetMap()) {
      a += e.second;
    }

    if (var.Int["verbose_stages"]) {
      std::cout << std::fixed;
      auto& mp = mt_.GetMap();
      std::cout << "mem=" << (sysinfo::GetMem() / double(1 << 20)) << " MB"
                << std::endl;
      ParseReport(mp, std::cout);
    }

    const auto& m = mk.begin()->second->GetMesh();
    const size_t nc = m.GetGlobalSize().prod();
    const size_t nt = var.Int["max_step"];
    const size_t ni = var.Int("iter", 1);

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
    std::cout << std::setprecision(5) << std::scientific;
    std::cout << "cells = " << nc << "\n"
              << "steps = " << nt << "\n"
              << "iters = " << ni << "\n"
              << "total = " << int(a) << " s"
              << " = ";
    std::cout << std::setfill('0');
    std::cout << std::setw(2) << h[0] << ":" << std::setw(2) << h[1] << ":"
              << std::setw(2) << h[2] << "." << std::setw(3) << h[3] << "\n";
    std::cout << "time/cell/iter = " << a / (nc * ni) << " s" << std::endl;
  }
}

template <class M>
void DistrMesh<M>::ReportOpenmp() {
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
          std::cout << "thread=" << std::setw(2) << omp_get_thread_num();
          std::cout << std::setw(8) << " cpu=" << std::setw(2)
                    << sched_getcpu();
          std::cout << std::endl;
        }
      }
    }
  }
#endif
}
