#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>
#include <utility>
#include <tuple>
#include <stdexcept>

#include "CubismDistr/Vars.h"
#include "CubismDistr/Interp.h"
#include "CubismDistr/Kernel.h"
#include "CubismDistr/KernelMesh.h"
//#include "CubismDistr/ICubism.h"
#include "CubismDistr/ILocal.h"
#include "CubismDistr/Distr.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/linear.hpp"
#include "hydro/advection.hpp"

#include "hydro/output.hpp"

#include "cmp.h"


template <class M>
class Advection : public KernelMesh<M> {
 public:
  using KM = KernelMesh<M>;
  using Mesh = M;
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using IdxCell = geom::IdxCell;
  static constexpr size_t dim = M::dim;

  Advection(Vars& par, const MyBlockInfo& bi, 
            std::function<Scal(Vect)> fu0,
            std::function<Vect(Vect)> fvel);
  void Run() override;

 protected:
  using KM::par;
  using KM::m;
  using KM::IsRoot;

 private:
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldCell<Scal> fc_src_;
  using AS = solver::AdvectionSolverExplicit<M>;
  std::unique_ptr<AS> as_; // advection solver
};

template <class _M>
class AdvectionFactory : public KernelMeshFactory<_M> {
 public:
  using M = _M;
  using K = Advection<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  AdvectionFactory(std::function<Scal(Vect)> fu0,
                   std::function<Vect(Vect)> fvel) 
      : fu0_(fu0), fvel_(fvel) {}
  K* Make(Vars& par, const MyBlockInfo& bi) override {
    return new K(par, bi, fu0_, fvel_);
  }

 private:
  std::function<Scal(Vect)> fu0_;
  std::function<Vect(Vect)> fvel_;
};

template <class M>
Advection<M>::Advection(Vars& par, const MyBlockInfo& bi,
                        std::function<Scal(Vect)> fu0,
                        std::function<Vect(Vect)> fvel)
    : KernelMesh<M>(par, bi)
{
  // initial field for advection
  geom::FieldCell<Scal> u0(m, 0);
  for (auto c : m.AllCells()) {
    Vect x = m.GetCenter(c);
    u0[c] = fu0(x);
  }

  // boundary conditions for advection (empty)
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> bc;

  // cell conditions for advection (empty)
  geom::MapCell<std::shared_ptr<solver::ConditionCell>> mc_cond;

  // flux
  ff_flux_.Reinit(m, 0);
  for (auto f : m.AllFaces()) {
    Vect x = m.GetCenter(f);
    ff_flux_[f] = fvel(x).dot(m.GetSurface(f));
  }
  
  // source
  fc_src_.Reinit(m, 0.);

  as_.reset(new AS(m, u0, bc, &ff_flux_, 
            &fc_src_, 0., par.Double["dt"], 
            0., 0., 1.));
}

#if 0
template <class M>
void Advection<M>::TestSolve(
    std::function<Scal(Vect)> fi /*initial*/,
    std::function<Scal(Vect)> fe /*exact*/,
    size_t cg /*check gap (separate from boundary)*/,
    std::string name,
    bool check /*abort if different from exact*/) {
  auto sem = m.GetSem("TestSolve");
  auto& bc = m.GetBlockCells();
  for (size_t n = 0; n < par.Int["num_steps"]; ++n) {
    if (sem.Nested("as->StartStep()")) {
      as_->StartStep();
    }
    if (sem.Nested("as->MakeIteration")) {
      as_->MakeIteration();
    }
    if (sem.Nested("as->FinishStep()")) {
      as_->FinishStep();
    }
  }
  if (sem("check")) {
    geom::GBlockCells<dim> cbc(MIdx(cg), gs_ - MIdx(2 * cg)); // check block
    FieldCell<bool> mask(m, false);
    for (auto i : m.AllCells()) {
      if (cbc.IsInside(bc.GetMIdx(i))) {
        mask[i] = true;
      }
    }
    auto& fc = as_->GetField();
    PFCMP(Mean(fc_exact_, m, mask), Mean(fc, m, mask), check);
    PFCMP(DiffMax(fc_exact_, fc, m, mask), 0., check);
  }
  if (sem("comm")) {
    fc_ = as_->GetField();
    m.Comm(&fc_);
    m.Comm(&fc_exact_);
  }
}
#endif

template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("advection");
  if (sem.Nested("start")) {
    as_->StartStep();
  }
  if (sem.Nested("iter")) {
    as_->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as_->FinishStep();
  }
  if (sem("comm")) { // comm for GetGlobalField
    auto& u = const_cast<geom::FieldCell<Scal>&>(as_->GetField());
    m.Comm(&u);
  }
}

#if 0
  par.Double["extent"] = 1.; // TODO don't overwrite extent
  Scal extent = par.Double["extent"];
  {
    MIdx p(par.Int["px"], par.Int["py"], par.Int["pz"]);
    MIdx b(par.Int["bx"], par.Int["by"], par.Int["bz"]);
    MIdx bs(par.Int["bsx"], par.Int["bsy"], par.Int["bsz"]);
    gs_ = p * b * bs;
  }
  Scal dx = extent / gs_.norminf(); 
  assert(dx > 0. && dx < extent);
  ge_ = Vect(gs_) * dx;
  Scal cfl = par.Double["cfl"];
  Vect vel(par.Vect["vel"]);
  Scal ns = par.Int["num_steps"];
  Scal nc = cfl * ns; // distance in cells passed
  Scal dt = dx * cfl / vel.norminf();
  par.Double.Set("dt", dt);
  Scal t = dt * ns;

  auto sem = m.GetSem("Run");
  auto f = [](Vect x) { 
      return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]); 
    };
  auto f0 = [](Vect x) { return 0.; };
  auto f1 = [](Vect x) { return 1.; };
  auto fx = [](Vect x) { return Vect(1., -1., 1.).dot(x); };
  auto fex = [=](Vect x) { x -= vel * t; return fx(x); };

  bool fatal = par.Int["fatal"];

  if (sem.Nested()) {
    TestSolve(f0, f0, 0, "f0", fatal);
  }
  if (sem.Nested()) {
    TestSolve(f1, f1, 0, "f1", fatal);
  }
  if (sem.Nested()) {
    TestSolve(fx, fex, 4, "fx", fatal);
  }
  if (sem.Nested()) {
    TestSolve(f, f, 0, "f", false);
  }
}

using Scal = double;
using M = geom::MeshStructured<Scal, 3>;
using K = Advection<M>;
using D = DistrMesh<KernelMeshFactory<M>>;
using FC = geom::FieldCell<Scal>;
using IdxCell = geom::IdxCell;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

#endif

// Dump values on inner cells to text file. 
// Format:
// <Nx> <Ny> <Nz>
// <data:x=0,y=0,z=0> <data:x=1,y=0,z=0> ...
template <class Scal, class B=geom::GBlock<geom::IdxCell, 3>>
void Dump(const geom::FieldCell<Scal>& u, B& b, std::string on) {
  std::ofstream o(on.c_str());

  auto s = b.GetDimensions();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : b.Range()) {
    o << u[c] << " ";
  }

  o << std::endl;
}

template <class M_>
class DistrSolver : public solver::UnsteadySolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using AS = solver::AdvectionSolverExplicit<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using IdxCell = geom::IdxCell;

  using KF = AdvectionFactory<M>;
  using D = DistrMesh<KernelMeshFactory<M>>; 

  DistrSolver(MPI_Comm comm, Vars apar, 
              std::function<Scal(Vect)> fu0,
              std::function<Vect(Vect)> fvel)
      : UnsteadySolver(0., apar.Double["dt"])
      , par(apar)
  {
    // Create kernel factory
    KF kf(fu0, fvel);

    // Initialize blocks
    Distr* b; // base
    if (par.Int["loc"]) {
      b = CreateLocal(comm, kf, par);
    } else {
      //b = CreateCubism(comm, kf, par);
    }

    d_.reset(dynamic_cast<D*>(b));
    if (!d_) {
      throw std::runtime_error("DistrSolver: Can't cast to D");
    }
  }
  void StartStep() {
    par.Double.Set("dt", this->GetTimeStep()); 
    par.Double.Set("t", this->GetTime()); 
    d_->Run();
  }
  geom::GBlock<IdxCell, dim> GetBlock() const {
    return d_->GetGlobalBlock();
  }
  geom::FieldCell<Scal> GetField() const {
    return d_->GetGlobalField(0);
  }

 private:
  Vars par;
  std::unique_ptr<D> d_;
};

void Main(MPI_Comm comm, Vars& par) {
  using M = geom::MeshStructured<double, 3>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto fu0 = [](Vect x) -> Scal { 
    return Vect(0.5,0.7,0.).dist(x) < 0.15 ? 1. : 0.; };
  auto fvelsincos = [](Vect x) -> Vect { 
    return Vect(-std::cos(x[1]) * std::sin(x[0]), 
                std::cos(x[0]) * std::sin(x[1]),
                0.); };
  auto fveluni = [](Vect x) -> Vect { 
    return Vect(-1., -1., 0.); };
  auto fvel = fveluni;

  DistrSolver<M> ds(comm, par, fu0, fvel);

  for (int k = 0; k < 3; ++k) {
    // single step for initial dump
    ds.StartStep();
    ds.FinishStep();

    auto b = ds.GetBlock();
    auto f = ds.GetField();
    Dump(f, b, "u" + std::to_string(k) + ".dat");

    for (int i = 0; i < 10; ++i) {
      ds.StartStep();
      ds.FinishStep();
    }
  }
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
      Main(comm, par);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, par);
  }

  MPI_Finalize();	
  return 0;
}
