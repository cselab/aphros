#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>

#include "CubismDistr/Vars.h"
#include "CubismDistr/Interp.h"
#include "CubismDistr/Kernel.h"
#include "CubismDistr/KernelMesh.h"
#include "CubismDistr/ICubism.h"
#include "CubismDistr/ILocal.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/linear.hpp"


template <class M>
class Simple : public KernelMesh<M> {
 public:
  using KM = KernelMesh<M>;
  using Mesh = M;
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using IdxCell = geom::IdxCell;
  static constexpr size_t dim = M::dim;

  Simple(Vars& par, const MyBlockInfo& bi);
  void Run() override;

 protected:
  using KM::par;
  using KM::bi_;
  using KM::m;

 private:
  void TestComm();
  void TestSolve();

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
class SimpleFactory : public KernelMeshFactory<_M> {
 public:
  using M = _M;
  using K = Simple<M>;
  K* Make(Vars& par, const MyBlockInfo& bi) override {
    return new K(par, bi);
  }
};

template <class M>
Simple<M>::Simple(Vars& par, const MyBlockInfo& bi) 
  : KernelMesh<M>(par, bi)
{}

bool Cmp(double a, double b) {
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
  MIdx bs(par.Int["bsx"], par.Int["bsy"], par.Int["bsz"]);
  {
    MIdx p(par.Int["px"], par.Int["py"], par.Int["pz"]);
    MIdx b(par.Int["bx"], par.Int["by"], par.Int["bz"]);
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
  if (sem.Nested("TestComm")) {
    TestComm();
  }
  if (sem.Nested("TestSolve")) {
    TestSolve();
  }
}


void Main(MPI_Comm comm, bool loc, Vars& par) {
  // read config files, parse arguments, maybe init global fields
  using Scal = double;
  using M = geom::MeshStructured<Scal, 3>;
  using K = Simple<M>;
  using KF = SimpleFactory<M>;
  using D = Distr;
  
  KF kf;

  // Initialize buffer mesh and make Simple for each block.
  std::unique_ptr<Distr> d;
  if (loc) {
    d.reset(CreateLocal(comm, kf, par));
  } else {
    d.reset(CreateCubism(comm, kf, par));
  }

  d->Run();
}

int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool isroot = (!rank);

  Vars par;   // parameter storage
  Interp ip(par); // interpretor (parser)

  std::string fn = "a.conf";
  if (argc == 1) {
    // nop
  } else if (argc == 2) {
    fn = argv[1];
  } else {
    if (isroot) {
      std::cerr << "usage: " << argv[0] << " [a.conf]" << std::endl;
    }
    return 1;
  }

  if (isroot) {
    std::cerr << "Loading config from '" << fn << "'" << std::endl;
  }

  std::ifstream f(fn);  // config file
  // Read file and run all commands
  ip.RunAll(f);   

  // Print vars on root
  if (isroot) {            
    std::cerr << "\n=== config begin ===" << std::endl;
    ip.PrintAll(std::cerr);
    std::cerr << "=== config end ===\n" << std::endl;
  }

  bool loc = par.Int["loc"];

  MPI_Comm comm;
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
