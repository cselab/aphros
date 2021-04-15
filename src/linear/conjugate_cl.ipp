// Created by Petr Karnakov on 03.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "linear.h"
#include "conjugate_cl.h"
#include "opencl/opencl.h"
#include "util/format.h"

DECLARE_FORCE_LINK_TARGET(linear_conjugate_cl);

namespace linear {

template <class M>
struct SolverConjugateCL<M>::Imp {
  using Owner = SolverConjugateCL<M>;

  Imp(Owner* owner, const Extra& extra_, const M&)
      : owner_(owner), conf(owner_->conf), extra(extra_) {}
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    auto sem = m.GetSem(__func__);
    struct Shared {
      FieldCell<Scal> fcu;
      FieldCell<Scal> fcr;
      FieldCell<Scal> fcp;
      FieldCell<Scal> fclp; // linear fc_system operator applied to p
      Scal dot_p_lp;
      Scal dot_r;
      Scal dot_r_prev;
      Scal max_r;
      FieldCell<Expr> fc_system;
    };
    struct {
      std::unique_ptr<Shared> shared_obj;
      Shared* shared;
      Scal dummy0;
      Scal dummy1;
      int iter = 0;
      Info info;
    } * ctx(sem);
    auto& t = *ctx;
    auto& ms = m.GetShared();
    if (sem("shared")) {
      if (m.IsLead()) {
        t.shared_obj = std::make_unique<Shared>();
        t.shared = t.shared_obj.get();
      }
      m.BcastFromLead(&t.shared);
    }
    if (sem("init0")) {
      auto& s = *t.shared;
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
      auto& s = *t.shared;
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
      auto& s = *t.shared;
      s.fcp = s.fcr;
      s.fclp.Reinit(ms, 0);
    }
    sem.LoopBegin();
    if (sem("iter0") && m.IsLead()) {
      auto& s = *t.shared;
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
      auto& s = *t.shared;
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
      auto& s = *t.shared;
      for (auto c : ms.Cells()) {
        s.fcp[c] = s.fcr[c] + (s.dot_r / (s.dot_r_prev + 1e-100)) * s.fcp[c];
      }
      ms.Comm(&s.fcp, M::CommStencil::direct_one);
    }
    if (sem("check")) {
      auto& s = *t.shared;
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
      auto& s = *t.shared;
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
  Owner* owner_;
  Conf& conf;
  Extra extra;
};

template <class M>
SolverConjugateCL<M>::SolverConjugateCL(
    const Conf& conf_, const Extra& extra, const M& m)
    : Base(conf_), imp(new Imp(this, extra, m)) {}

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
        this->GetConf(var, prefix), extra, m);
  }
};

} // namespace linear
