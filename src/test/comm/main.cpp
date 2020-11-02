// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <mpi.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "distr/distrsolver.h"
#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "linear/linear.h"
#include "solver/pois.h"
#include "solver/solver.h"
#include "util/linear.h"
#include "util/suspender.h"

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

struct GPar {};

template <class M_>
class Simple : public KernelMeshPar<M_, GPar> {
 public:
  using P = KernelMeshPar<M_, GPar>; // parent
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = GPar;
  static constexpr size_t dim = M::dim;

  using P::P;
  void Run() override;

 protected:
  using P::bi_;
  using P::m;
  using P::var;

 private:
  void TestComm();
  void TestReduce();
  void TestScatter();
  void TestPois();

  // TODO: revise 1e-10
  bool Cmp(Scal a, Scal b) {
    return std::abs(a - b) < 1e-10;
  }
  bool Cmp(Vect a, Vect b) {
    return a.dist(b) < 1e-10;
  }
  bool Cmp(MIdx a, MIdx b) {
    return a == b;
  }
  template <class T>
  bool Cmp(const std::vector<T>& a, const std::vector<T>& b) {
    return a == b;
  }

  FieldCell<Scal> fc_;
  FieldCell<Vect> fcv_;
  Scal r_; // test Reduce
  std::pair<Scal, int> rsi_; // test Reduce minloc
  std::vector<Scal> rvs_; // reduction vector<Scal> (concatenation)
  std::vector<int> rvi_; // reduction vector<int> (concatenation)
  std::vector<std::vector<int>> rvvi_; // reduction vector<vector<int>>
  std::vector<std::vector<Scal>> rvvs_; // reduction vector<vector<int>>
};

template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r = std::max(r, std::abs(u[i] - v[i]));
  }
  return r;
}

template <class Idx, class M>
typename M::Scal DiffMean(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r += std::abs(u[i] - v[i]);
    w += 1.;
  }
  return r / w;
}

template <class Idx, class M>
typename M::Scal Max(const GField<typename M::Scal, Idx>& u, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r = std::max(r, u[i]);
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Mean(const GField<typename M::Scal, Idx>& u, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0.;
  for (auto i : m.template GetRangeIn<Idx>()) {
    r += u[i];
    w += 1.;
  }
  return r / w;
}

#define CMP(a, b)                                                \
  if (!Cmp(a, b)) {                                              \
    throw std::runtime_error(                                    \
        std::string(__FILE__) + ":" + std::to_string(__LINE__)); \
  }

// Print CMP
#define PCMP(a, b)                                                        \
  {                                                                       \
    auto _fmt = std::cerr.flags();                                        \
    std::cerr << std::scientific << std::setprecision(16) << "CMP " << #a \
              << "=" << a << "\n    " << #b << "=" << b << std::endl;     \
    CMP(a, b);                                                            \
    std::cerr.flags(_fmt);                                                \
  }

// Print CMP standard flags
#define PCMPF(a, b)                                                  \
  std::cerr << "CMP " << std::setw(15) << #a << "=" << a << "\n    " \
            << std::setw(15) << #b << "=" << b << std::endl;         \
  CMP(a, b);

template <class M>
void Simple<M>::TestComm() {
  auto sem = m.GetSem("Comm");
  auto f = [](Vect v) {
    for (size_t i = 0; i < dim; ++i) {
      while (v[i] < 0.) {
        v[i] += 1.;
      }
      while (v[i] > 1.) {
        v[i] -= 1.;
      }
    }
    return std::sin(v[0]) * std::cos(v[1]) * std::exp(v[2]);
  };
  auto fv = [=](Vect v) { return Vect(f(v), f(v * 2.), f(v * 3.)); };
  auto& bc = m.GetIndexCells();
  if (sem("init")) {
    fc_.Reinit(m);
    fcv_.Reinit(m);
    for (auto c : m.Cells()) {
      auto x = m.GetCenter(c);
      fc_[c] = f(x);
      fcv_[c] = fv(x);
    }
    m.Comm(&fc_);
    m.Comm(&fcv_);
  }
  if (sem("check")) {
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      if (!Cmp(fc_[c], f(x))) {
        std::cerr << bc.GetMIdx(c) << " " << fc_[c] << " != " << f(x) << " "
                  << std::endl;
        throw std::runtime_error("");
      }
      if (!Cmp(fcv_[c], fv(x))) {
        std::cerr << bc.GetMIdx(c) << " " << fcv_[c] << " != " << fv(x) << " "
                  << std::endl;
        throw std::runtime_error("");
      }
    }
  }
}

template <class M>
void Simple<M>::TestReduce() {
  auto sem = m.GetSem("Reduce");
  MIdx proc(var.Int["px"], var.Int["py"], var.Int["pz"]);
  MIdx block(var.Int["bx"], var.Int["by"], var.Int["bz"]);
  GBlock<IdxCell, dim> bq(proc * block);
  GIndex<size_t, dim> indexc(proc * block);
  GBlock<IdxCell, dim> qp(proc);
  GBlock<IdxCell, dim> qb(block);
  auto f = [](MIdx w) {
    return std::sin(w[0] + 1.7) * std::cos(w[1]) * std::exp(w[2] * 0.1);
  };
  if (sem("sum")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, Reduction::sum);
  }
  if (sem("sum-check")) {
    Scal exact = 0.;
    for (auto b : bq) {
      exact += f(b);
    }
    PCMP(r_, exact);
  }
  if (sem("prod")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, Reduction::prod);
  }
  if (sem("prod-check")) {
    Scal exact = 1.;
    for (auto b : bq) {
      exact *= f(b);
    }
    PCMP(r_, exact);
  }
  if (sem("max")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, Reduction::max);
  }
  if (sem("max-check")) {
    Scal s = -std::numeric_limits<Scal>::max();
    for (auto b : bq) {
      s = std::max(s, f(b));
    }
    PCMP(r_, s);
  }
  if (sem("min")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, Reduction::min);
  }
  if (sem("min-check")) {
    Scal s = std::numeric_limits<Scal>::max();
    for (auto b : bq) {
      s = std::min(s, f(b));
    }
    PCMP(r_, s);
  }
  if (sem("minloc")) {
    MIdx w(bi_.index);
    rsi_ = std::make_pair(f(w), indexc.GetIdx(w));
    m.Reduce(&rsi_, Reduction::minloc);
  }
  if (sem("minloc-check")) {
    Scal s = std::numeric_limits<Scal>::max();
    MIdx ws;
    for (auto w : bq) {
      if (f(w) < s) {
        ws = w;
        s = f(w);
      }
    }
    PCMP(rsi_.first, s);
    PCMP(indexc.GetMIdx(rsi_.second), ws);
  }
  if (sem("maxloc")) {
    MIdx w(bi_.index);
    rsi_ = std::make_pair(f(w), indexc.GetIdx(w));
    m.Reduce(&rsi_, Reduction::maxloc);
  }
  if (sem("maxloc-check")) {
    Scal s = -std::numeric_limits<Scal>::max();
    MIdx ws;
    for (auto w : bq) {
      if (f(w) > s) {
        ws = w;
        s = f(w);
      }
    }
    PCMP(rsi_.first, s);
    PCMP(indexc.GetMIdx(rsi_.second), ws);
  }
  if (sem("cat")) {
    MIdx w(bi_.index);
    rvs_.resize(0);
    rvs_.push_back(indexc.GetIdx(w));
    m.Reduce(&rvs_, Reduction::concat);
  }
  if (sem("cat-check")) {
    std::vector<Scal> exact;
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = block * wp + wb;
        exact.push_back(indexc.GetIdx(w));
      }
    }
    if (m.IsRoot()) {
      std::cout << "****" << NAMEVALUE(m.GetId()) << std::endl;
      std::sort(rvs_.begin(), rvs_.end());
      std::sort(exact.begin(), exact.end());
      PCMPF(rvs_, exact);
    }
  }
  const size_t q = 10;
  if (sem("cati")) {
    MIdx w(bi_.index);
    rvi_.resize(0);
    size_t i = indexc.GetIdx(w);
    for (size_t j = 0; j < i % q; ++j) {
      rvi_.push_back(q * i + j);
    }
    m.Reduce(&rvi_, Reduction::concat);
  }
  if (sem("cati-check")) {
    std::vector<int> rvi;
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = block * wp + wb;
        size_t i = indexc.GetIdx(w);
        for (size_t j = 0; j < i % q; ++j) {
          rvi.push_back(q * i + j);
        }
      }
    }
    if (m.IsRoot()) {
      std::sort(rvi_.begin(), rvi_.end());
      std::sort(rvi.begin(), rvi.end());
      PCMP(rvi_, rvi);
    }
  }
  if (sem("catvi")) {
    MIdx w(bi_.index);
    rvvi_.resize(0);
    size_t i = indexc.GetIdx(w);
    for (size_t j = 0; j < i % q; ++j) {
      rvvi_.push_back(std::vector<int>({int(q * i + j)}));
    }
    m.Reduce(&rvvi_, Reduction::concat);
  }
  if (sem("catvi-check")) {
    std::vector<std::vector<int>> rvvi;
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = block * wp + wb;
        size_t i = indexc.GetIdx(w);
        for (size_t j = 0; j < i % q; ++j) {
          rvvi.push_back(std::vector<int>({int(q * i + j)}));
        }
      }
    }
    if (m.IsRoot()) {
      std::sort(rvvi_.begin(), rvvi_.end());
      std::sort(rvvi.begin(), rvvi.end());
      PCMP(rvvi_, rvvi);
    }
  }
  if (sem("bcast-catvi")) {
    MIdx w(bi_.index);
    rvvi_.resize(0);
    size_t i = indexc.GetIdx(w);
    for (size_t j = 0; j < (i + 5) % q; ++j) {
      rvvi_.push_back(std::vector<int>({int(100 * i + j)}));
    }
    m.Bcast(&rvvi_);
  }
  if (sem("bcast-catvi-check")) {
    // init for root block
    std::vector<std::vector<int>> rvvi;
    size_t i = indexc.GetIdx(MIdx(0));
    for (size_t j = 0; j < (i + 5) % q; ++j) {
      rvvi.push_back(std::vector<int>({int(100 * i + j)}));
    }
    PCMP(rvvi_, rvvi);
  }
}

template <class M>
void Simple<M>::TestScatter() {
  auto sem = m.GetSem("scatter");
  MIdx proc(var.Int["px"], var.Int["py"], var.Int["pz"]);
  MIdx block(var.Int["bx"], var.Int["by"], var.Int["bz"]);
  GBlock<IdxCell, dim> bq(proc * block);
  GIndex<size_t, dim> indexc(proc * block);
  GBlock<IdxCell, dim> qp(proc);
  GBlock<IdxCell, dim> qb(block);
  auto GetBlockData = [](size_t i) {
    std::vector<Scal> r;
    r.push_back(Scal(i));
    for (size_t j = 0; j < i; ++j) {
      r.push_back(Scal(j));
    }
    return r;
  };
  if (sem("vector<Scal>")) {
    rvvs_.clear();
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = block * wp + wb;
        size_t i = indexc.GetIdx(w);
        rvvs_.push_back(GetBlockData(i));
      }
    }
    if (m.IsRoot()) {
      std::cerr << "rvvs=" << rvvs_ << std::endl;
      m.Scatter({&rvvs_, &rvs_});
    } else {
      m.Scatter({nullptr, &rvs_});
    }
  }
  if (sem("check")) {
    assert(rvs_.size() > 0);
    PCMPF(rvs_, GetBlockData(std::lround(rvs_[0])));
  }
  if (sem("gather")) {
    MIdx w(bi_.index);
    size_t i = indexc.GetIdx(w);
    rvvs_ = {GetBlockData(i)};
    m.Reduce(&rvvs_, Reduction::concat);
  }
  if (sem("scatter")) {
    if (m.IsRoot()) {
      std::cerr << "rvvs=" << rvvs_ << std::endl;
      m.Scatter({&rvvs_, &rvs_});
    } else {
      m.Scatter({nullptr, &rvs_});
    }
  }
  if (sem("check")) {
    MIdx w(bi_.index);
    size_t i = indexc.GetIdx(w);
    PCMPF(rvs_, GetBlockData(i));
  }
}

template <class M>
void Simple<M>::TestPois() {
  auto sem = m.GetSem("pois");
  struct {
    MapEmbed<BCond<Scal>> mebc;
    FieldCell<Scal> fc_sol;
    FieldCell<Scal> fc_rhs;
    FieldCell<Scal> fc_exsol;
    std::shared_ptr<linear::Solver<M>> linsolver;
  } * ctx(sem);
  auto& t = *ctx;
  // exact solution
  auto fe = [](Vect x) {
    Scal pi = M_PI;
    Scal k = 2 * pi;
    return std::sin(x[0] * k) * std::sin(x[1] * k) * std::sin(x[2] * k);
  };
  // rhs
  auto fr = [=](Vect x, Vect h) {
    auto f = fe;
    Vect h0(h[0], 0., 0.);
    Vect h1(0., h[1], 0.);
    Vect h2(0., 0., h[2]);
    Scal d0 = (f(x - h0) - 2. * f(x) + f(x + h0)) / sqr(h[0]);
    Scal d1 = (f(x - h1) - 2. * f(x) + f(x + h1)) / sqr(h[1]);
    Scal d2 = (f(x - h2) - 2. * f(x) + f(x + h2)) / sqr(h[2]);
    return d0 + d1 + d2;
  };

  if (sem("init")) {
    if (m.IsRoot()) {
      std::cout << "\n\n*** TestPois() ***" << std::endl;
    }
    t.linsolver = ULinear<M>::MakeLinearSolver(var, "vort");
    // exact solution
    t.fc_exsol.Reinit(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      t.fc_exsol[i] = fe(x);
    }
    // rhs
    t.fc_rhs.Reinit(m);
    for (auto c : m.AllCells()) {
      t.fc_rhs[c] = fr(m.GetCenter(c), m.GetCellSize());
    }
  }
  if (sem.Nested("solve")) {
    SolvePoisson(t.fc_sol, t.fc_rhs, t.mebc, t.linsolver, m);
  }
  if (sem("check")) {
    // Check
    PCMP(Mean(t.fc_exsol, m), Mean(t.fc_sol, m));
    PCMP(DiffMean(t.fc_exsol, t.fc_sol, m), 0.);
    PCMP(DiffMax(t.fc_exsol, t.fc_sol, m), 0.);
  }
}

template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem("Run");
  if (sem.Nested()) {
    TestComm();
  }
  if (sem.Nested()) {
    TestReduce();
  }
  if (sem.Nested()) {
    TestScatter();
  }
  if (sem.Nested()) {
    TestPois();
  }
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using K = Simple<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
