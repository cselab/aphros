// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "debug/isnan.h"
#include "debug/linear.h"
#include "geom/mesh.h"


// V is expression: v[0] * c + v[1] * cxm + ... v[6] * czp + v[7]
// lx: initial guess
template <
    class M, class Scal = typename M::Scal,
    class V = generic::Vect<Scal, M::dim * 2 + 2>>
typename M::LS ConvertLsCompact(
    const FieldCell<V>& fce, std::vector<Scal>& la, std::vector<Scal>& lb,
    std::vector<Scal>& lx, const M& m) {
  using LS = typename M::LS;
  LS l;

  // stencil
  // XXX: adhoc, assumed order in m.GetCell(c)
  /*
  l.st.emplace_back(0, 0, 0);   // 0
  l.st.emplace_back(-1, 0, 0);  // 1
  l.st.emplace_back(1, 0, 0);   // 2
  l.st.emplace_back(0, -1, 0);  // 3
  l.st.emplace_back(0, 1, 0);   // 4
  l.st.emplace_back(0, 0, -1);  // 5
  l.st.emplace_back(0, 0, 1);   // 6
  */

  l.st.emplace_back(0, 0, -1); // 5
  l.st.emplace_back(0, -1, 0); // 3
  l.st.emplace_back(-1, 0, 0); // 1
  l.st.emplace_back(0, 0, 0); // 0
  l.st.emplace_back(1, 0, 0); // 2
  l.st.emplace_back(0, 1, 0); // 4
  l.st.emplace_back(0, 0, 1); // 6
  assert(l.st.size() == V::dim - 1);

  size_t n = m.GetInBlockCells().size(); // number of equations
  la.resize(n * l.st.size());
  lb.resize(n);
  lx.resize(n);

  // fill matrix coeffs
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      la[i++] = e[5];
      la[i++] = e[3];
      la[i++] = e[1];
      la[i++] = e[0];
      la[i++] = e[2];
      la[i++] = e[4];
      la[i++] = e[6];
    }
    assert(i == n * l.st.size());
  }

  // fill rhs and zero solution
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      lb[i] = -e[V::dim - 1];
      ++i;
    }
    assert(i == lb.size());
  }

  l.a = &la;
  l.b = &lb;
  l.x = &lx;
  return l;
}

// Solve linear system fc_system = 0
// fc_system: expressions [i]
// fc_init: initial guess
// prefix: custom refix for solver config
// Output:
// fc_sol: solution [a]
// m.GetSolveTmp(): modified temporary fields
template <class M>
void Solve(
    const FieldCell<typename M::Expr>& fc_system,
    const FieldCell<typename M::Scal>* fc_init,
    FieldCell<typename M::Scal>& fc_sol, typename M::LS::T type, M& m,
    std::string prefix = "") {
  using Scal = typename M::Scal;
  auto sem = m.GetSem("solve");
  if (type == M::LS::T::symm && m.flags.check_symmetry) {
    if (sem.Nested()) {
      UDebug<M>::CheckSymmetry(fc_system, m);
    }
  }
  if (sem("solve")) {
    std::vector<Scal>* lsa;
    std::vector<Scal>* lsb;
    std::vector<Scal>* lsx;
    m.GetSolveTmp(lsa, lsb, lsx);
    lsx->resize(m.GetInBlockCells().size());
    if (fc_init) {
      size_t i = 0;
      for (auto c : m.Cells()) {
        (*lsx)[i++] = (*fc_init)[c];
      }
    } else {
      size_t i = 0;
      for (auto c : m.Cells()) {
        (void)c;
        (*lsx)[i++] = 0;
      }
    }
    auto l = ConvertLsCompact(fc_system, *lsa, *lsb, *lsx, m);
    l.t = type;
    l.prefix = prefix;
    m.Solve(l);
  }
  if (sem("copy")) {
    std::vector<Scal>* lsa;
    std::vector<Scal>* lsb;
    std::vector<Scal>* lsx;
    m.GetSolveTmp(lsa, lsb, lsx);

    fc_sol.Reinit(m);
    size_t i = 0;
    for (auto c : m.Cells()) {
      fc_sol[c] = (*lsx)[i++];
    }
    m.Comm(&fc_sol);

    if (m.flags.linreport && m.IsRoot()) {
      std::cout << std::scientific;
      std::cout << "linear solver '" + fc_system.GetName() + "':"
                << " res=" << m.GetResidual() << " iter=" << m.GetIter()
                << std::endl;
    }
  }
}


namespace linear {

template <class M>
class Solver {
 public:
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;

  struct Conf {
    Scal tol = 0;
    int maxiter = 100;
  };

  struct Info {
    Scal residual;
    int iter;
  };

  Solver(const Conf& conf_) : conf(conf_) {}
  virtual ~Solver() = default;
  // Solves linear system
  //   system(x) = 0
  //
  // Input:
  // fc_system: coefficients `e` of linear expressions
  //   e[0] * x[c] + sum_i(e[i + 1] * x[c(i)]) + e.back() = 0
  // fc_init: initial guess, may equal &fc_sol, nullptr to assume zero guess.
  //
  // Output:
  // fc_sol: solution
  // Info: final residual and number of iterations
  virtual Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) = 0;
  virtual void SetConf(const Conf& c) {
    conf = c;
  }
  virtual const Conf& GetConf() {
    return conf;
  }

 protected:
  Conf conf;
};

template <class M>
class SolverHypre : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra { // extra config
    std::string solver = "pcg"; // name of the solver to use
    int print = 0; // print level, 0 for none
  };
  SolverHypre(const Conf& conf, const Extra& extra);
  ~SolverHypre();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
class SolverConjugate : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {};
  SolverConjugate(const Conf& conf, const Extra& extra);
  ~SolverConjugate();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
class SolverJacobi : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {};
  SolverJacobi(const Conf& conf, const Extra& extra);
  ~SolverJacobi();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
class SolverDefault : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {
    typename M::LS::T type;
    std::string prefix = "";
  };
  SolverDefault(const Conf& conf, const Extra& extra);
  ~SolverDefault();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

} // namespace linear

