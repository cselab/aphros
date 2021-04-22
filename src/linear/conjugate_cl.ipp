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
      auto& cl = s.cl;
      s.d_system.EnqueueWrite(cl.queue, s.fc_system.data()->data());
      s.d_fcu.EnqueueWrite(cl.queue, s.fcu.data());

      s.kernel_linear.EnqueueWithArgs(
          cl.queue, cl.global_size, cl.local_size, cl.start, cl.lead_y,
          cl.lead_z, s.d_system, s.d_fcu, s.d_fcr, Scal(-1), Scal(-1));
    }
    if (sem.Nested()) {
      s.cl.Comm(m, s.d_fcr);
    }
    if (sem("init2") && m.IsLead()) {
      auto& cl = s.cl;
      s.d_fcp.EnqueueCopyFrom(cl.queue, s.d_fcr);
    }
    sem.LoopBegin();
    if (sem("iter0_cl") && m.IsLead()) {
      auto& cl = s.cl;
      s.kernel_linear.EnqueueWithArgs(
          cl.queue, cl.global_size, cl.local_size, cl.start, cl.lead_y,
          cl.lead_z, s.d_system, s.d_fcp, s.d_fclp, Scal(1), Scal(0));
    }
    if (sem("iter0") && m.IsLead()) {
      auto& cl = s.cl;
      s.dot_r_prev = cl.Dot(s.d_fcr, s.d_fcr);
      s.dot_p_lp = cl.Dot(s.d_fcp, s.d_fclp);

      ms.Reduce(&s.dot_r_prev, Reduction::sum);
      ms.Reduce(&s.dot_p_lp, Reduction::sum);
    }
    if (sem("iter1") && m.IsLead()) {
      auto& cl = s.cl;
      const Scal alpha = s.dot_r_prev / (s.dot_p_lp + 1e-100);
      s.kernel_accum.EnqueueWithArgs(
          cl.queue, cl.global_size, cl.local_size, cl.start, cl.lead_y,
          cl.lead_z, alpha, s.d_fcp, s.d_fcu);
      s.kernel_accum.EnqueueWithArgs(
          cl.queue, cl.global_size, cl.local_size, cl.start, cl.lead_y,
          cl.lead_z, -alpha, s.d_fclp, s.d_fcr);

      s.dot_r = cl.Dot(s.d_fcr, s.d_fcr);
      s.max_r = cl.Max(s.d_fcr);

      ms.Reduce(&s.dot_r, Reduction::sum);
      ms.Reduce(&s.max_r, Reduction::max);
    }
    if (sem("iter2_cl") && m.IsLead()) {
      auto& cl = s.cl;
      s.kernel_iter2.EnqueueWithArgs(
          cl.queue, cl.global_size, cl.local_size, cl.start, cl.lead_y,
          cl.lead_z, s.d_fcp, s.d_fcr, s.dot_r, s.dot_r_prev, s.d_fcp);
      cl.queue.Finish();
    }
    if (sem.Nested()) {
      s.cl.Comm(m, s.d_fcp);
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
    if (sem("result") && m.IsLead()) {
      auto& cl = s.cl;
      s.d_fcu.EnqueueRead(cl.queue, s.fcu.data());
      cl.queue.Finish();
    }
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
      d_fcu.Create(cl.context, cl.size);
      d_fcp.Create(cl.context, cl.size);
      d_fclp.Create(cl.context, cl.size);
      d_fcr.Create(cl.context, cl.size);
      d_system.Create(cl.context, cl.size * (M::dim * 2 + 2));
      program.CreateFromString(kProgram, cl.context, cl.device);
      kernel_iter2.Create(program, "iter2");
      kernel_linear.Create(program, "linear");
      kernel_accum.Create(program, "accum");
    }
    FieldCell<Scal> fcu;
    Scal dot_p_lp;
    Scal dot_r;
    Scal dot_r_prev;
    Scal max_r;
    FieldCell<Expr> fc_system;
    OpenCL<M> cl;
    using Kernel = typename OpenCL<M>::Kernel;
    using Program = typename OpenCL<M>::Program;
    template <class T>
    using Buffer = typename OpenCL<M>::template Buffer<T>;
    Program program;
    Kernel kernel_iter2;
    Kernel kernel_linear;
    Kernel kernel_accum;
    Buffer<Scal> d_fcu;
    Buffer<Scal> d_fcp;
    Buffer<Scal> d_fclp; // linear fc_system operator applied to p
    Buffer<Scal> d_fcr;
    Buffer<Scal> d_system;
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
