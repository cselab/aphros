#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>

#include "CubismDistr/Vars.h"
#include "CubismDistr/Interp.h"
#include "CubismDistr/Kernel.h"
#include "CubismDistr/Cubism.h"
#include "CubismDistr/Local.h"
#include "CubismDistr/Hydro.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/linear.hpp"


template <class M>
class Simple : public Kernel {
 public:
  using Mesh = M;
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using IdxCell = typename geom::IdxCell;
  static constexpr size_t dim = M::dim;

  Simple(Vars& par, const MyBlockInfo& bi);
  void Run() override;
  M& GetMesh() { return m; }

 private:
  void TestComm();
  void TestSolve();

  Vars& par;
  std::string name_;
  MyBlockInfo bi_;
  M m;
  geom::FieldCell<Scal> fc_;
  // LS
  using Expr = solver::Expression<Scal, IdxCell, 1 + dim * 2>;
  geom::FieldCell<Expr> fc_system_;
  std::vector<Scal> lsa_;
  std::vector<Scal> lsb_;
  std::vector<Scal> lsx_;
  geom::FieldCell<Scal> fc_sol_;
  geom::FieldCell<Scal> fc_exsol_;
};

template <class _M>
class SimpleFactory : public KernelFactory {
 public:
  using M = _M;
  using K = Simple<M>;
  std::unique_ptr<Kernel> Make(Vars& par, const MyBlockInfo& bi) override {
    return std::unique_ptr<K>(new K(par, bi));
  }
};

template <class M>
Simple<M>::Simple(Vars& par, const MyBlockInfo& bi) 
  : par(par), bi_(bi), m(CreateMesh<M>(bi))
{
  name_ = 
      "[" + std::to_string(bi.index[0]) +
      "," + std::to_string(bi.index[1]) +
      "," + std::to_string(bi.index[2]) + "]";
}

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-10;
}

template <class Idx, class M>
typename M::Scal DiffMax(
    const geom::GField<typename M::Scal, Idx>& u,
    const geom::GField<typename M::Scal, Idx>& v,
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
    const geom::GField<typename M::Scal, Idx>& u,
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
    const geom::GField<typename M::Scal, Idx>& u,
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
  assert(Cmp(a, b)); 

// Print CMP
#define PCMP(a, b) \
  std::cerr \
    << std::scientific << std::setprecision(16) \
    << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b); 


template <class M>
void Simple<M>::TestComm() {
  auto sem = m.GetSem("TestComm");
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
        assert(false);
      }
    }
  }
}

template <class M>
void Simple<M>::TestSolve() {
  auto sem = m.GetSem("TestSolve");
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
  {
    MIdx p(par.Int["px"], par.Int["py"], par.Int["pz"]);
    MIdx b(par.Int["bx"], par.Int["by"], par.Int["bz"]);
    using B = MyBlock;
    MIdx s(B::sx, B::sy, B::sz); // block size inner
    gs = p * b * s;
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

      bool per = par.Int["periodic"];


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

    using LS = typename Mesh::LS;
    LS l;
    // Get stencil from first inner cell
    {
      IdxCell c = *m.Cells().begin(); 
      auto& e = fc_system_[c];
      for (size_t j = 0; j < e.size(); ++j) {
        MIdx dm = bc.GetMIdx(e[j].idx) - bc.GetMIdx(c);
        l.st.emplace_back(dm);
      }
    }

    int bs = _BLOCKSIZE_;
    int n = bs * bs *bs;
    using MIdx = typename Mesh::MIdx;
    lsa_.resize(n*l.st.size());
    lsb_.resize(n, 1.);
    lsx_.resize(n, 0.);

    // fill matrix coeffs
    {
      size_t i = 0;
      for (auto c : m.Cells()) {
        auto& e = fc_system_[c];
        for (size_t j = 0; j < e.size(); ++j) {
          // Check stencil
          if (e[j].idx != bc.GetIdx(bc.GetMIdx(c) + MIdx(l.st[j]))) {
            std::cerr << "***"
                << " MIdx(c)=" << bc.GetMIdx(c)
                << " MIdx(e[j].idx)=" << bc.GetMIdx(e[j].idx)
                << " l.st[j]=" << MIdx(l.st[j]) 
                << std::endl;
            assert(false);
          }
          lsa_[i] = e[j].coeff;
          ++i;
        }
      }
      assert(i == n * l.st.size());
    }

    // fill rhs and zero solution
    {
      size_t i = 0;
      for (auto c : m.Cells()) {
        auto& e = fc_system_[c];
        lsb_[i] = -e.GetConstant();
        lsx_[i] = 0.;
        ++i;
      }
      assert(i == lsb_.size());
    }

    l.a = &lsa_;
    l.b = &lsb_;
    l.x = &lsx_;
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
  if (sem.Nested("TestComm")) {
    TestComm();
  }
  if (sem.Nested("TestSolve")) {
    TestSolve();
  }
}


void Main(MPI_Comm comm, bool loc, Vars& par) {
  // read config files, parse arguments, maybe init global fields
  using M = geom::MeshStructured<Scal, 3>;
  using K = Simple<M>;
  using KF = SimpleFactory<M>;
  using D = Distr;
  
  KF kf;

  const int es = 8;
  const int hl = par.Int["hl"];
  const int bs = 16;
  
  // Initialize buffer mesh and make Simple for each block.
  std::unique_ptr<Distr> d;
  if (loc) {
    d.reset(new Local<KF>(comm, kf, bs, es, hl, par));
  } else {
    d.reset(new Cubism<KF>(comm, kf, bs, es, hl, par));
  }

  while (!d->IsDone()) {
    d->Step();
  }
}


int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm comm;

  Vars par;
  Interp ip(par);

  std::ifstream f("a.conf");
  ip.RunAll(f);
  if (rank == 0) {
    ip.PrintAll();
  }

  bool loc = par.Int["loc"];

  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      Main(comm, loc, par);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, loc, par);
  }


  MPI_Finalize();	
  return 0;
}
