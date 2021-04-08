// Created by Petr Karnakov on 03.10.2020
// Copyright 2020 ETH Zurich

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "hypre.h"
#include "linear_hypre.h"

DECLARE_FORCE_LINK_TARGET(linear_hypre);

namespace linear {

template <class M>
struct SolverHypre<M>::Imp {
  using Owner = SolverHypre<M>;

  Imp(Owner* owner, const Extra& extra_, const M&)
      : owner_(owner), conf(owner_->conf), extra(extra_) {}
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    using MIdx = typename M::MIdx;
    auto sem = m.GetSem(__func__);
    struct {
      // local on all blocks
      std::vector<MIdx> stencil;
      std::vector<Scal> data_a;
      std::vector<Scal> data_b;
      std::vector<Scal> data_x;
      Info info;

      // reduced to lead block
      std::vector<MIdx> origin;
      std::vector<MIdx> size;
      std::vector<std::vector<Scal>*> ptr_a;
      std::vector<std::vector<Scal>*> ptr_b;
      std::vector<std::vector<Scal>*> ptr_x;
      std::vector<Info*> ptr_info;
      std::unique_ptr<Hypre> hypre;
    } * ctx(sem);
    auto& t = *ctx;
    if (m.flags.check_symmetry &&
        (extra.solver == "pcg" || extra.solver == "pcg+smg")) {
      if (sem.Nested()) {
        UDebug<M>::CheckSymmetry(fc_system, m);
      }
    }
    if (sem("solve")) {
      const auto bic = m.GetInBlockCells();
      t.origin.push_back(bic.GetBegin());
      t.size.push_back(bic.GetSize());

      // copy data from l to block-local buffer
      t.stencil = {
          MIdx{0, 0, 0}, MIdx{-1, 0, 0}, MIdx{1, 0, 0}, MIdx{0, -1, 0},
          MIdx{0, 1, 0}, MIdx{0, 0, -1}, MIdx{0, 0, 1},
      };

      t.data_a.resize(t.stencil.size() * bic.size());
      t.data_b.resize(bic.size());
      t.data_x.resize(bic.size(), 0);

      { // matrix coeffs
        size_t i = 0;
        for (auto c : m.Cells()) {
          for (size_t k = 0; k < 7; ++k) {
            t.data_a[i++] = fc_system[c][k];
          }
        }
      }

      { // rhs
        size_t i = 0;
        for (auto c : m.Cells()) {
          t.data_b[i++] = -fc_system[c].back();
        }
      }

      if (fc_init) { // initial guess
        size_t i = 0;
        for (auto c : m.Cells()) {
          t.data_x[i++] = (*fc_init)[c];
        }
      }

      // pass pointers to block-local data to the lead block
      t.ptr_a.push_back(&t.data_a);
      t.ptr_b.push_back(&t.data_b);
      t.ptr_x.push_back(&t.data_x);
      t.ptr_info.push_back(&t.info);

      m.GatherToLead(&t.origin);
      m.GatherToLead(&t.size);
      m.GatherToLead(&t.ptr_a);
      m.GatherToLead(&t.ptr_b);
      m.GatherToLead(&t.ptr_x);
      m.GatherToLead(&t.ptr_info);
    }
    if (sem("hypre") && m.IsLead()) {
      using HypreBlock = typename Hypre::Block;
      using HypreMIdx = typename Hypre::MIdx;

      const size_t nblocks = t.origin.size();
      fassert_equal(t.size.size(), nblocks);
      fassert_equal(t.ptr_a.size(), nblocks);
      fassert_equal(t.ptr_b.size(), nblocks);
      fassert_equal(t.ptr_x.size(), nblocks);

      std::vector<HypreBlock> blocks(nblocks);
      for (size_t i = 0; i < nblocks; ++i) {
        HypreBlock& block = blocks[i];
        block.l = t.origin[i];
        block.u = t.origin[i] + t.size[i] - MIdx(1);
        for (auto w : t.stencil) {
          block.stencil.push_back(w);
        }
        block.a = t.ptr_a[i];
        block.r = t.ptr_b[i];
        block.x = t.ptr_x[i];
      }

      const HypreMIdx per = MIdx(m.flags.is_periodic);
      t.hypre = std::make_unique<Hypre>(
          m.GetMpiComm(), blocks, m.GetGlobalSize(), per);

      t.hypre->Solve(conf.tol, extra.print, extra.solver, conf.maxiter);

      // fill local info and update by pointers on other blocks
      t.info.residual = t.hypre->GetResidual();
      t.info.iter = t.hypre->GetIter();
      for (size_t i = 0; i < nblocks; ++i) {
        t.ptr_info[i] = &t.info;
      }
    }
    if (sem()) {
      // copy solution from t.data_x to field
      fc_sol.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc_sol[c] = t.data_x[i++];
      }
      m.Comm(&fc_sol);
      if (m.flags.linreport && m.IsRoot()) {
        std::cerr << std::scientific;
        std::cerr << "linear(hypre) '" + fc_system.GetName() + "':"
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
SolverHypre<M>::SolverHypre(const Conf& conf_, const Extra& extra, const M& m)
    : Base(conf_), imp(new Imp(this, extra, m)) {}

template <class M>
SolverHypre<M>::~SolverHypre() = default;

template <class M>
auto SolverHypre<M>::Solve(
    const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, M& m) -> Info {
  return imp->Solve(fc_system, fc_init, fc_sol, m);
}

template <class M>
class ModuleLinearHypre : public ModuleLinear<M> {
 public:
  ModuleLinearHypre() : ModuleLinear<M>("hypre") {}
  std::unique_ptr<Solver<M>> Make(
      const Vars& var, std::string prefix, const M& m) override {
    auto addprefix = [prefix](std::string name) {
      return "hypre_" + prefix + "_" + name;
    };
    typename SolverHypre<M>::Extra extra;
    extra.solver = var.String[addprefix("solver")];
    extra.print = var.Int["hypre_print"];
    return std::make_unique<linear::SolverHypre<M>>(
        this->GetConf(var, prefix), extra, m);
  }
};

using M = MeshCartesian<double, 3>;

bool kReg_hypre[] = {
    RegisterModule<ModuleLinearHypre<M>>(),
};

} // namespace linear
