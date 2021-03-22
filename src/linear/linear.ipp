// Created by Petr Karnakov on 03.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "linear.h"

DECLARE_FORCE_LINK_TARGET(linear_conjugate);
DECLARE_FORCE_LINK_TARGET(linear_jacobi);

namespace linear {

template <class M>
struct SolverConjugate<M>::Imp {
  using Owner = SolverConjugate<M>;

  Imp(Owner* owner, const Extra& extra_, const M&)
      : owner_(owner), conf(owner_->conf), extra(extra_) {}
  // TODO: use fc_init
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      FieldCell<Scal> fcu;
      FieldCell<Scal> fcr;
      FieldCell<Scal> fcp;
      FieldCell<Scal> fc_opp; // linear fc_system operator applied to p
      Scal dot_p_opp;
      Scal dot_r;
      Scal dot_r_prev;
      Scal max_r;

      int iter = 0;
      Info info;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("init")) {
      if (fc_init) {
        t.fcu = *fc_init;
      } else {
        t.fcu.Reinit(m, 0);
      }
      t.fcr.Reinit(m);
      for (auto c : m.Cells()) {
        const auto& e = fc_system[c];
        Scal u = t.fcu[c] * e[0] + e.back();
        for (auto q : m.Nci(c)) {
          u += t.fcu[m.GetCell(c, q)] * e[1 + q.raw()];
        }
        t.fcr[c] = -u;
      }
      m.Comm(&t.fcr, M::CommStencil::direct_one);
    }
    if (sem("init")) {
      t.fcp = t.fcr;
      t.fc_opp.Reinit(m);
    }
    sem.LoopBegin();
    if (sem("iter")) {
      for (auto c : m.Cells()) {
        const auto& e = fc_system[c];
        Scal p = t.fcp[c] * e[0];
        for (auto q : m.Nci(c)) {
          p += t.fcp[m.GetCell(c, q)] * e[1 + q.raw()];
        }
        t.fc_opp[c] = p;
      }

      t.dot_r_prev = 0;
      t.dot_p_opp = 0;
      for (auto c : m.Cells()) {
        t.dot_r_prev += sqr(t.fcr[c]);
        t.dot_p_opp += t.fcp[c] * t.fc_opp[c];
      }
      m.Reduce(&t.dot_r_prev, Reduction::sum);
      m.Reduce(&t.dot_p_opp, Reduction::sum);
    }
    if (sem("iter2")) {
      const Scal alpha = t.dot_r_prev / (t.dot_p_opp + 1e-100);
      t.dot_r = 0;
      t.max_r = 0;
      for (auto c : m.Cells()) {
        t.fcu[c] += alpha * t.fcp[c];
        t.fcr[c] -= alpha * t.fc_opp[c];
        t.dot_r += sqr(t.fcr[c]);
        t.max_r = std::max(t.max_r, std::abs(t.fcr[c]));
      }
      m.Reduce(&t.dot_r, Reduction::sum);
      m.Reduce(&t.max_r, Reduction::max);
    }
    if (sem("iter3")) {
      for (auto c : m.Cells()) {
        t.fcp[c] = t.fcr[c] + (t.dot_r / (t.dot_r_prev + 1e-100)) * t.fcp[c];
      }
      m.Comm(&t.fcp, M::CommStencil::direct_one);
    }
    if (sem("check")) {
      if (extra.residual_max) {
        t.info.residual = t.max_r / m.GetCellSize().prod();
      } else { // L2-norm
        t.info.residual = std::sqrt(t.dot_r / m.GetCellSize().prod());
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
      fc_sol = t.fcu;
      m.Comm(&fc_sol, M::CommStencil::direct_one);
      if (m.flags.linreport && m.IsRoot()) {
        std::cerr << std::scientific;
        std::cerr << "linear(conjugate) '" + fc_system.GetName() + "':"
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
SolverConjugate<M>::SolverConjugate(
    const Conf& conf_, const Extra& extra, const M& m)
    : Base(conf_), imp(new Imp(this, extra, m)) {}

template <class M>
SolverConjugate<M>::~SolverConjugate() = default;

template <class M>
auto SolverConjugate<M>::Solve(
    const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, M& m) -> Info {
  return imp->Solve(fc_system, fc_init, fc_sol, m);
}

template <class M>
struct SolverJacobi<M>::Imp {
  using Owner = SolverJacobi<M>;

  Imp(Owner* owner, const Extra& extra_, const M&)
      : owner_(owner), conf(owner_->conf), extra(extra_) {}
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      FieldCell<Scal> fcu;
      FieldCell<Scal> fcu_new;
      Scal maxdiff;
      int iter = 0;
      Info info;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("init")) {
      if (fc_init) {
        t.fcu = *fc_init;
      } else {
        t.fcu.Reinit(m, 0);
      }
      t.fcu_new.Reinit(m, 0);
    }
    sem.LoopBegin();
    if (sem("iter")) {
      t.maxdiff = 0;
      for (auto c : m.Cells()) {
        const auto& e = fc_system[c];
        Scal nondiag = e.back();
        for (auto q : m.Nci(c)) {
          nondiag += t.fcu[m.GetCell(c, q)] * e[1 + q.raw()];
        }
        t.fcu_new[c] = -nondiag / e[0];
        t.maxdiff = std::max(t.maxdiff, std::abs(t.fcu_new[c] - t.fcu[c]));
      }
      t.fcu.swap(t.fcu_new);
      m.Comm(&t.fcu);
      m.Reduce(&t.maxdiff, Reduction::max);
    }
    if (sem("check")) {
      t.info.residual = t.maxdiff;
      ++t.iter;
      t.info.iter = t.iter;
      if (t.iter >= conf.miniter &&
          (t.iter > conf.maxiter || t.info.residual < conf.tol)) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (sem("result")) {
      fc_sol = t.fcu;
      m.Comm(&fc_sol);
      if (m.flags.linreport && m.IsRoot()) {
        std::cerr << std::scientific;
        std::cerr << "linear(jacobi) '" + fc_system.GetName() + "':"
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
SolverJacobi<M>::SolverJacobi(const Conf& conf_, const Extra& extra, const M& m)
    : Base(conf_), imp(new Imp(this, extra, m)) {}

template <class M>
SolverJacobi<M>::~SolverJacobi() = default;

template <class M>
auto SolverJacobi<M>::Solve(
    const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, M& m) -> Info {
  return imp->Solve(fc_system, fc_init, fc_sol, m);
}

template <class M>
class ModuleLinearConjugate : public ModuleLinear<M> {
 public:
  ModuleLinearConjugate() : ModuleLinear<M>("conjugate") {}
  std::unique_ptr<Solver<M>> Make(
      const Vars& var, std::string prefix, const M& m) override {
    auto addprefix = [prefix](std::string name) {
      return "linsolver_" + prefix + "_" + name;
    };
    typename linear::SolverConjugate<M>::Extra extra;
    extra.residual_max = var.Int[addprefix("maxnorm")];
    return std::make_unique<linear::SolverConjugate<M>>(
        this->GetConf(var, prefix), extra, m);
  }
};

template <class M>
class ModuleLinearJacobi : public ModuleLinear<M> {
 public:
  ModuleLinearJacobi() : ModuleLinear<M>("jacobi") {}
  std::unique_ptr<Solver<M>> Make(
      const Vars& var, std::string prefix, const M& m) override {
    typename linear::SolverJacobi<M>::Extra extra;
    return std::make_unique<linear::SolverJacobi<M>>(
        this->GetConf(var, prefix), extra, m);
  }
};

} // namespace linear
