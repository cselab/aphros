#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>

#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
#include "util/suspender.h"
#include "solver/solver.h"
#include "linear/linear.h"

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

  FieldCell<Scal> fc_;
  // LS
  using Expr = solver::Expression<Scal, IdxCell, 1 + dim * 2>;
  FieldCell<Expr> fc_system_;
  std::vector<Scal> lsa_;
  std::vector<Scal> lsb_;
  std::vector<Scal> lsx_;
  FieldCell<Scal> fc_sol_;
  FieldCell<Scal> fc_exsol_;
  Scal r_; // test Reduce
};

bool Cmp(double a, double b) {
  return std::abs(a - b) < 1e-10;
}

template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    r = std::max(r, std::abs(u[i] - v[i]));
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Max(
    const GField<typename M::Scal, Idx>& u,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
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
  for (auto i : m.template Get<Idx>()) {
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
    for (int i = 0; i < dim; ++i) {
      while (v[i] < 0.) {
        v[i] += 1.;
      }
      while (v[i] > 1.) {
        v[i] -= 1.;
      }
    }
    return std::sin(v[0]) * std::cos(v[1]) * std::exp(v[2]); 
  };
  auto& bc = m.GetBlockCells();
  if (sem("init")) {
    fc_.Reinit(m);
    for (auto i : m.Cells()) {
      fc_[i] = f(m.GetCenter(i));
    }
    m.Comm(&fc_);
  }
  if (sem("check")) {
    for (auto i : m.AllCells()) {
      auto x = m.GetCenter(i);
      if (!Cmp(fc_[i], f(m.GetCenter(i)))) {
        std::cerr 
          << bc.GetMIdx(i) << " " 
          << fc_[i] << " != " << f(x) << " "
          << std::endl;
        throw std::runtime_error("");
      }
    }
  }
}

template <class M>
void Simple<M>::TestReduce() {
  auto sem = m.GetSem("Reduce");
  auto f = [](MIdx w) {
    return std::sin(w[0]) * std::cos(w[1]) * std::exp(w[2]); 
  };
  if (sem("init")) {
    r_ = f(MIdx(bi_.index));
    m.Reduce(&r_, "sum");
  }
  if (sem("check")) {
    MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
    MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
    GBlock<size_t, dim> bb(p * b);
    Scal s = 0.;
    for (auto b : bb) {
      s += f(b);
    }
    PCMP(r_, s);
  }
}


template <class M>
void Simple<M>::TestSolve() {
  auto sem = m.GetSem("Solve");
  auto& bc = m.GetBlockCells();
  auto f = [](Vect v) { 
    for (int i = 0; i < dim; ++i) {
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

      Vect x = m.GetCenter(i);
      Vect xpx = m.GetCenter(ipx);
      Vect xmx = m.GetCenter(imx);
      Vect xpy = m.GetCenter(ipy);
      Vect xmy = m.GetCenter(imy);
      Vect xpz = m.GetCenter(ipz);
      Vect xmz = m.GetCenter(imz);

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

      /*
      if (i == *m.Cells().begin()) {
        e *= 0.;
      }
      */

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
