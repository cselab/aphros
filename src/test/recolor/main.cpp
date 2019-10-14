#undef NDEBUG
#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <limits>
#include <cmath>

#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
#include "util/suspender.h"
#include "solver/solver.h"
#include "linear/linear.h"
#include "solver/pois.h"

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
  using P::var;
  using P::bi_;
  using P::m;

 private:
  void TestComm();
  void TestSolve();
  void TestReduce();
  void TestPois();

  // TODO: revise 1e-10
  bool Cmp(Scal a, Scal b) { return std::abs(a - b) < 1e-10; }
  bool Cmp(Vect a, Vect b) { return a.dist(b) < 1e-10; }
  bool Cmp(MIdx a, MIdx b) { return a == b; }
  template <class T>
  bool Cmp(const std::vector<T>& a, 
           const std::vector<T>& b) { return a == b; }

  FieldCell<Scal> fc_;
  FieldCell<Vect> fcv_;
  // LS
  using Expr = solver::Expression<Scal, IdxCell, 1 + dim * 2>;
  FieldCell<Expr> fc_system_;
  std::vector<Scal> lsa_;
  std::vector<Scal> lsb_;
  std::vector<Scal> lsx_;
  FieldCell<Scal> fc_sol_;
  FieldCell<Scal> fc_exsol_;
  FieldCell<Scal> fcr_; // rhs
  Scal r_; // test Reduce
  std::pair<Scal, int> rsi_; // test Reduce minloc
  std::vector<Scal> rvs_; // reduction vector<Scal> (concatenation)
  std::vector<int> rvi_; // reduction vector<int> (concatenation)
  std::vector<std::vector<int>> rvvi_; // reduction vector<vector<int>> 
  std::shared_ptr<solver::PoisSolver<M>> ps_;
};

template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetIn<Idx>()) {
    r = std::max(r, std::abs(u[i] - v[i]));
  }
  return r;
}

template <class Idx, class M>
typename M::Scal DiffMean(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0;
  for (auto i : m.template GetIn<Idx>()) {
    r += std::abs(u[i] - v[i]);
    w += 1.;
  }
  return r / w;
}

template <class Idx, class M>
typename M::Scal Max(
    const GField<typename M::Scal, Idx>& u,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetIn<Idx>()) {
    r = std::max(r, u[i]);
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Mean(
    const GField<typename M::Scal, Idx>& u,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0.;
  for (auto i : m.template GetIn<Idx>()) {
    r += u[i];
    w += 1.;
  }
  return r / w;
}


#define CMP(a, b) \
  if (!Cmp(a,b)) { \
    throw std::runtime_error(std::string(__FILE__) + std::to_string(__LINE__)); \
  }

// Print CMP
#define PCMP(a, b) \
  std::cerr \
    << std::scientific << std::setprecision(16) \
    << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
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
  auto fv = [=](Vect v) { 
    return Vect(f(v), f(v * 2.), f(v * 3.));
  };
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
        std::cerr 
          << bc.GetMIdx(c) << " " 
          << fc_[c] << " != " << f(x) << " "
          << std::endl;
        throw std::runtime_error("");
      }
      if (!Cmp(fcv_[c], fv(x))) {
        std::cerr 
          << bc.GetMIdx(c) << " " 
          << fcv_[c] << " != " << fv(x) << " "
          << std::endl;
        throw std::runtime_error("");
      }
    }
  }
}

template <class M>
void Simple<M>::TestReduce() {
  auto sem = m.GetSem("Reduce");
  MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
  MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
  GBlock<IdxCell, dim> bq(p * b); 
  GIndex<size_t, dim> ndq(p * b); 
  GBlock<IdxCell, dim> qp(p); 
  GBlock<IdxCell, dim> qb(b); 
  auto f = [](MIdx w) {
    return std::sin(w[0]+1.7) * std::cos(w[1]) * std::exp(w[2]); 
  };
  if (sem("sum")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, "sum");
  }
  if (sem("sum-check")) {
    Scal s = 0.;
    for (auto b : bq) {
      s += f(b);
    }
    PCMP(r_, s);
  }
  if (sem("prod")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, "prod");
  }
  if (sem("prod-check")) {
    Scal s = 1.;
    for (auto b : bq) {
      s *= f(b);
    }
    PCMP(r_, s);
  }
  if (sem("max")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, "max");
  }
  if (sem("max-check")) {
    Scal s = std::numeric_limits<Scal>::min();
    for (auto b : bq) {
      s = std::max(s, f(b));
    }
    PCMP(r_, s);
  }
  if (sem("min")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, "min");
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
    rsi_ = std::make_pair(f(w), ndq.GetIdx(w));
    m.Reduce(std::make_shared<typename M::OpMinloc>(&rsi_));
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
    PCMP(ndq.GetMIdx(rsi_.second), ws);
  }
  if (sem("maxloc")) {
    MIdx w(bi_.index);
    rsi_ = std::make_pair(f(w), ndq.GetIdx(w));
    m.Reduce(std::make_shared<typename M::OpMaxloc>(&rsi_));
  }
  if (sem("maxloc-check")) {
    Scal s = std::numeric_limits<Scal>::min();
    MIdx ws;
    for (auto w : bq) {
      if (f(w) > s) {
        ws = w;
        s = f(w);
      }
    }
    PCMP(rsi_.first, s);
    PCMP(ndq.GetMIdx(rsi_.second), ws);
  }
  if (sem("cat")) {
    MIdx w(bi_.index);
    rvs_.resize(0);
    rvs_.push_back(ndq.GetIdx(w));
    using T = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<T>(&rvs_));
  }
  if (sem("cat-check")) {
    std::vector<Scal> rvs;
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = b * wp + wb;
        rvs.push_back(ndq.GetIdx(w));
      }
    }
    if (m.IsRoot()) {
      PCMP(rvs_, rvs);
    }
  }
  const size_t q = 10;
  if (sem("cati")) {
    MIdx w(bi_.index);
    rvi_.resize(0);
    size_t i = ndq.GetIdx(w);
    for (size_t j = 0; j < i % q; ++j) {
      rvi_.push_back(q * i + j);
    }
    using T = typename M::template OpCatT<int>;
    m.Reduce(std::make_shared<T>(&rvi_));
  }
  if (sem("cati-check")) {
    std::vector<int> rvi;
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = b * wp + wb;
        size_t i = ndq.GetIdx(w);
        for (size_t j = 0; j < i % q; ++j) {
          rvi.push_back(q * i + j);
        }
      }
    }
    if (m.IsRoot()) {
      PCMP(rvi_, rvi);
    }
  }
  if (sem("catvi")) {
    MIdx w(bi_.index);
    rvvi_.resize(0);
    size_t i = ndq.GetIdx(w);
    for (size_t j = 0; j < i % q; ++j) {
      rvvi_.push_back(std::vector<int>({int(q * i + j)}));
    }
    using T = typename M::template OpCatVT<int>;
    m.Reduce(std::make_shared<T>(&rvvi_));
  }
  if (sem("catvi-check")) {
    std::vector<std::vector<int>> rvvi;
    for (auto wp : qp) {
      for (auto wb : qb) {
        auto w = b * wp + wb;
        size_t i = ndq.GetIdx(w);
        for (size_t j = 0; j < i % q; ++j) {
          rvvi.push_back(std::vector<int>({int(q * i + j)}));
        }
      }
    }
    if (m.IsRoot()) {
      PCMP(rvvi_, rvvi);
    }
  }
  if (sem("bcast-catvi")) {
    MIdx w(bi_.index);
    rvvi_.resize(0);
    size_t i = ndq.GetIdx(w);
    for (size_t j = 0; j < (i + 5) % q; ++j) {
      rvvi_.push_back(std::vector<int>({int(100 * i + j)}));
    }
    using T = typename M::template OpCatVT<int>;
    m.Bcast(std::make_shared<T>(&rvvi_));
  }
  if (sem("bcast-catvi-check")) {
    // init for root block
    std::vector<std::vector<int>> rvvi;
    size_t i = ndq.GetIdx(MIdx(0));
    for (size_t j = 0; j < (i + 5) % q; ++j) {
      rvvi.push_back(std::vector<int>({int(100 * i + j)}));
    }
    PCMP(rvvi_, rvvi);
  }
}


template <class M>
void Simple<M>::TestSolve() {
  auto sem = m.GetSem("Solve");
  auto& bc = m.GetIndexCells();
  auto f = [](Vect v) { 
    for (size_t i = 0; i < dim; ++i) {
      while (v[i] < 0.) {
        v[i] += 1.;
      }
      while (v[i] > 1.) {
        v[i] -= 1.;
      }
    }
    return std::sin(v[0]* v[1]) * 
        std::cos(v[1] * v[2]) * 
        std::exp(v[2] + v[0]); 
  };

  // global mesh size
  MIdx gs;
  MIdx bs(var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]);
  {
    MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
    MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
    gs = p * b * bs;
  }

  if (sem("init")) {
    if (m.IsRoot()) {
      std::cout << "\n\n*** TestSolve() ***" << std::endl;
    }
    // exact solution
    fc_exsol_.Reinit(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_exsol_[i] = f(x);
    }
    // init hydro system
    fc_system_.Reinit(m);
    for (auto i : m.Cells()) {
      MIdx mi = bc.GetMIdx(i);
      IdxCell ipx = bc.GetIdx(mi + MIdx(1, 0, 0));
      IdxCell imx = bc.GetIdx(mi + MIdx(-1, 0, 0));
      IdxCell ipy = bc.GetIdx(mi + MIdx(0, 1, 0));
      IdxCell imy = bc.GetIdx(mi + MIdx(0, -1, 0));
      IdxCell ipz = bc.GetIdx(mi + MIdx(0, 0, 1));
      IdxCell imz = bc.GetIdx(mi + MIdx(0, 0, -1));

      MIdx mpx = bc.GetMIdx(ipx);
      MIdx mmx = bc.GetMIdx(imx);
      MIdx mpy = bc.GetMIdx(ipy);
      MIdx mmy = bc.GetMIdx(imy);
      MIdx mpz = bc.GetMIdx(ipz);
      MIdx mmz = bc.GetMIdx(imz);

      auto& e = fc_system_[i];
      e.Clear();

      bool per = var.Int["periodic"];

      e.InsertTerm(mpx < gs || per       ? -0.0 : 0., ipx);
      e.InsertTerm(MIdx(0) <= mmx || per ? -0.0 : 0., imx);
      e.InsertTerm(mpy < gs || per       ? -1.0 : 0., ipy);
      e.InsertTerm(MIdx(0) <= mmy || per ? -1.0 : 0., imy);
      e.InsertTerm(mpz < gs || per       ? -1.0 : 0., ipz);
      e.InsertTerm(MIdx(0) <= mmz || per ? -1.0 : 0., imz);

      e.InsertTerm(6., i);

      Scal r = e.Evaluate(fc_exsol_);
      e.SetConstant(-r);
    }

    using LS = typename M::LS;
    LS l = ConvertLs(fc_system_, lsa_, lsb_, lsx_, m);
    m.Solve(l);
  }
  if (sem("check")) {
    // Copy solution to field
    fc_sol_.Reinit(m);
    size_t i = 0;
    for (auto c : m.Cells()) {
      fc_sol_[c] = lsx_[i++];
    }
    // Check
    PCMP(Mean(fc_exsol_, m), Mean(fc_sol_, m));
    PCMP(DiffMax(fc_exsol_, fc_sol_, m), 0.);
  }
}

template <class M>
void Simple<M>::TestPois() {
  auto sem = m.GetSem("pois");
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
    // solver
    MapFace<std::shared_ptr<solver::CondFace>> mf;
    ps_ = std::make_shared<solver::PoisSolver<M>>(mf, m);
    // exact solution
    fc_exsol_.Reinit(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_exsol_[i] = fe(x);
    }
    // rhs
    fcr_.Reinit(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fcr_[i] = fr(x, m.GetCellSize());
    }
  }
  if (sem.Nested("solve")) {
    ps_->Solve(fcr_);
  }
  if (sem("check")) {
    fc_sol_ = ps_->GetField();
    // Check
    PCMP(Mean(fc_exsol_, m), Mean(fc_sol_, m));
    PCMP(DiffMean(fc_exsol_, fc_sol_, m), 0.);
    PCMP(DiffMax(fc_exsol_, fc_sol_, m), 0.);
  }
}


template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem("Run");
  if (sem.Nested("Comm")) {
    TestComm();
  }
  if (sem.Nested("Solve")) {
    TestSolve();
  }
  if (sem.Nested("Reduce")) {
    TestReduce();
  }
  if (sem.Nested("Pois")) {
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
