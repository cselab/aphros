// Created by Petr Karnakov on 03.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "conjugate_cl.h"
#include "linear.h"
#include "opencl/opencl.h"
#include "parse/vars.h"
#include "util/format.h"

DECLARE_FORCE_LINK_TARGET(linear_conjugate_cl);

const char* kProgram =
#include "conjugate_cl.inc"
    ;

namespace linear {

template <class M>
struct SolverConjugateCL<M>::Imp {
  using Owner = SolverConjugateCL<M>;

  Imp(Owner* owner, const Extra& extra_, const M& m, const Vars& var)
      : owner_(owner), conf(owner_->conf), extra(extra_) {
    if (m.IsLead()) {
      shared_obj_ = std::make_unique<Shared>(m, var);
      shared = shared_obj_.get();
    }
  }
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      int iter = 0;
      Info info;
    } * ctx(sem);
    auto& t = *ctx;
    auto& ms = m.GetShared();
    if (sem()) {
      m.BcastFromLead(&shared);
    }
    auto& s = *shared;
    if (sem("init0")) {
      if (fc_init) {
        LocalToShared(*fc_init, s.fcu, m);
        if (m.IsLead()) {
          ms.Comm(&s.fcu, M::CommStencil::direct_one);
        }
      } else {
        if (m.IsLead()) {
          s.fcu.Reinit(ms, 0);
        }
      }
      LocalToShared(fc_system, s.fc_system, m);
    }
    if (sem("init1") && m.IsLead()) {
      s.fcr.Reinit(ms);
      for (auto c : ms.Cells()) {
        const auto& e = s.fc_system[c];
        Scal u = s.fcu[c] * e[0] + e.back();
        for (auto q : ms.Nci(c)) {
          u += s.fcu[ms.GetCell(c, q)] * e[1 + q.raw()];
        }
        s.fcr[c] = -u;
      }
      ms.Comm(&s.fcr, M::CommStencil::direct_one);
    }
    if (sem("init2") && m.IsLead()) {
      s.fcp = s.fcr;
      s.fclp.Reinit(ms, 0);
    }
    sem.LoopBegin();
    if (sem("iter0") && m.IsLead()) {
      for (auto c : ms.Cells()) {
        const auto& e = s.fc_system[c];
        Scal p = s.fcp[c] * e[0];
        for (auto q : ms.Nci(c)) {
          p += s.fcp[ms.GetCell(c, q)] * e[1 + q.raw()];
        }
        s.fclp[c] = p;
      }

      s.dot_r_prev = 0;
      s.dot_p_lp = 0;
      for (auto c : ms.Cells()) {
        s.dot_r_prev += sqr(s.fcr[c]);
        s.dot_p_lp += s.fcp[c] * s.fclp[c];
      }

      ms.Reduce(&s.dot_r_prev, Reduction::sum);
      ms.Reduce(&s.dot_p_lp, Reduction::sum);
    }
    if (sem("iter1") && m.IsLead()) {
      const Scal alpha = s.dot_r_prev / (s.dot_p_lp + 1e-100);
      s.dot_r = 0;
      s.max_r = 0;

      for (auto c : ms.Cells()) {
        s.fcu[c] += alpha * s.fcp[c];
        s.fcr[c] -= alpha * s.fclp[c];
        s.dot_r += sqr(s.fcr[c]);
        s.max_r = std::max(s.max_r, std::abs(s.fcr[c]));
      }

      ms.Reduce(&s.dot_r, Reduction::sum);
      ms.Reduce(&s.max_r, Reduction::max);
    }
    if (sem("iter2") && m.IsLead()) {
      for (auto c : ms.Cells()) {
        s.fcp[c] = s.fcr[c] + (s.dot_r / (s.dot_r_prev + 1e-100)) * s.fcp[c];
        s.fcp[c] += 1; // XXX
      }
      ms.Comm(&s.fcp, M::CommStencil::direct_one);
    }
    if (sem("iter2_cl") && m.IsLead()) {
      auto& cl = s.cl;
      auto info = OpenCL<M>::Device::GetDeviceInfo(cl.device.platform);
      s.d_fcp.EnqueueWrite(cl.queue, s.fcp.data());
      s.kernel.EnqueueWithArgs(
          cl.queue, cl.global_size, cl.local_size, cl.start, cl.lead_y,
          cl.lead_z, s.d_fcp, s.d_fcp);
      s.d_fcp.EnqueueRead(cl.queue, s.fcp.data());
      cl.queue.Finish();
    }
    if (sem("check")) {
      if (extra.residual_max) {
        t.info.residual = s.max_r / ms.GetCellSize().prod();
      } else { // L2-norm
        t.info.residual = std::sqrt(s.dot_r / ms.GetCellSize().prod());
      }
      ++t.iter;
      t.info.iter = t.iter;
      if (t.iter >= conf.miniter &&
          (t.iter > conf.maxiter || t.info.residual < conf.tol)) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (sem("result")) {
      SharedToLocal(s.fcu, fc_sol, m);
      m.Comm(&fc_sol, M::CommStencil::direct_one);
      if (m.flags.linreport && m.IsRoot()) {
        std::cerr << std::scientific;
        std::cerr << "linear(conjugate_cl) '" + fc_system.GetName() + "':"
                  << " res=" << t.info.residual << " iter=" << t.info.iter
                  << std::endl;
      }
    }
    if (sem()) {
    }
    return t.info;
  }

 private:
  struct Shared {
    Shared(const M& m, const Vars& var) : cl(m.GetShared(), var) {
      d_fcp.Create(cl.context, cl.size);
      program.CreateFromString(kProgram, cl.context, cl.device);
      kernel.Create(program, "iter2");
      CLCALL(clGetKernelWorkGroupInfo(
          kernel, cl.device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
          &max_work_size, NULL));
    }
    FieldCell<Scal> fcu;
    FieldCell<Scal> fcr;
    FieldCell<Scal> fcp;
    FieldCell<Scal> fclp; // linear fc_system operator applied to p
    Scal dot_p_lp;
    Scal dot_r;
    Scal dot_r_prev;
    Scal max_r;
    FieldCell<Expr> fc_system;
    OpenCL<M> cl;
    using Kernel = typename OpenCL<M>::Kernel;
    using Program = typename OpenCL<M>::Program;
    template <class T>
    using Buffer = typename OpenCL<M>::Buffer<T>;
    Program program;
    Kernel kernel;
    size_t max_work_size;
    Buffer<Scal> d_fcp;
  };

  Owner* owner_;
  Conf& conf;
  Extra extra;

  std::unique_ptr<Shared> shared_obj_;
  Shared* shared;
};

template <class M>
SolverConjugateCL<M>::SolverConjugateCL(
    const Conf& conf_, const Extra& extra, const M& m, const Vars& var)
    : Base(conf_), imp(new Imp(this, extra, m, var)) {}

template <class M>
SolverConjugateCL<M>::~SolverConjugateCL() = default;

template <class M>
auto SolverConjugateCL<M>::Solve(
    const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, M& m) -> Info {
  return imp->Solve(fc_system, fc_init, fc_sol, m);
}

template <class M>
class ModuleLinearConjugateCL : public ModuleLinear<M> {
 public:
  ModuleLinearConjugateCL() : ModuleLinear<M>("conjugate_cl") {}
  std::unique_ptr<Solver<M>> Make(
      const Vars& var, std::string prefix, const M& m) override {
    auto addprefix = [prefix](std::string name) {
      return "linsolver_" + prefix + "_" + name;
    };
    typename linear::SolverConjugateCL<M>::Extra extra;
    extra.residual_max = var.Int[addprefix("maxnorm")];
    return std::make_unique<linear::SolverConjugateCL<M>>(
        this->GetConf(var, prefix), extra, m, var);
  }
};

} // namespace linear
