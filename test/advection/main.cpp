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
#include "CubismDistr/KernelMeshPar.h"
//#include "CubismDistr/ICubism.h"
#include "CubismDistr/ILocal.h"
#include "CubismDistr/Distr.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/linear.hpp"
#include "hydro/advection.hpp"
#include "hydro/advection_vof.hpp"
#include "hydro/advection_tvd.hpp"
#include "hydro/init_vel.h"
#include "hydro/init_u.h"
#include "hydro/dump.h"

#include "hydro/output.hpp"

#include "cmp.h"

using namespace geom;

template <class M_>
struct GPar {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  std::function<Scal(Vect /*t*/)> fu0; // initial volume fraction
  std::function<Vect(Vect /*x*/,Scal /*t*/)> fvel; // velocity
};

template <class M_>
class Advection : public KernelMeshPar<M_, GPar<M_>> {
 public:
  using P = KernelMeshPar<M_, GPar<M_>>; // parent
  using M = M_;
  using Mesh = M;
  using Par = typename P::Par;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  //using P::P; // inherit constructor
  Advection(Vars& var, const MyBlockInfo& bi, Par& par);
  void Run() override;

 protected:
  using P::var;
  using P::par_;
  using P::m;
  using P::IsRoot;
  using P::IsLead;

 private:
  FieldFace<Scal> ff_flux_;
  FieldCell<Scal> fc_src_;
  using AS = solver::AdvectionSolver<M>;
  using AS1 = solver::AdvectionSolverExplicit<M>;
  using AS2 = solver::Vof<M>;
  std::unique_ptr<AS> as_; // advection solver
  std::function<Vect(Vect,Scal)> fvel_;
  Scal maxvel_; // maximum velocity relative to cell length [1/time]
                // cfl = dt * maxvel
  Scal sumu_; // sum of fluid volume
  typename AS1::Par aspar1_;
  typename AS2::Par aspar2_;
};

template <class M>
Advection<M>::Advection(Vars& var, const MyBlockInfo& bi, Par& par)
    : KernelMeshPar<M, Par>(var, bi, par)
{
  // initial field for advection
  FieldCell<Scal> u0(m, 0);
  for (auto c : m.AllCells()) {
    Vect x = m.GetCenter(c);
    u0[c] = par_.fu0(x);
  }

  // boundary conditions for advection (empty)
  MapFace<std::shared_ptr<solver::ConditionFace>> bc;

  // cell conditions for advection (empty)
  MapCell<std::shared_ptr<solver::ConditionCell>> mc_cond;

  // flux
  ff_flux_.Reinit(m, 0);
  for (auto f : m.AllFaces()) {
    Vect x = m.GetCenter(f);
    ff_flux_[f] = par_.fvel(x, 0.).dot(m.GetSurface(f));
  }
  
  // source
  fc_src_.Reinit(m, 0.);

  std::string as = var.String["advection_solver"];
  if (as == "tvd") {
    auto& p = aspar1_;
    p.sharp = var.Double["sharp"];
    p.sharpo = var.Double["sharpo"];
    p.split = var.Int["split"];
    as_.reset(new AS1(m, u0, bc, &ff_flux_, 
              &fc_src_, 0., var.Double["dt"], &p));
  } else if (as == "vof") {
    auto& p = aspar2_;
    as_.reset(new AS2(m, u0, bc, &ff_flux_, 
              &fc_src_, 0., var.Double["dt"], &p));
  } else {
    throw std::runtime_error("Unknown advection_solver=" + as);
  }
}


template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("advection");
  if (sem("comm")) { // comm for GetGlobalField, initial
    auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
    m.Comm(&u);
  }
  sem.LoopBegin();
  if (sem("empty")) {
    //nop // TODO: bugfix loop, empty stage needed
  }
  if (sem("checkloop")) {
    if (as_->GetTime() >= var.Double["tmax"]) {
      sem.LoopBreak();
    }
  }
  if (sem("vel")) {
    const Scal t = as_->GetTime();
    for (auto f : m.AllFaces()) {
      Vect x = m.GetCenter(f);
      ff_flux_[f] = par_.fvel(x, t).dot(m.GetSurface(f));
    }

    maxvel_ = 0.;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        auto f = m.GetNeighbourFace(c, q);
        maxvel_ = std::max(maxvel_, std::abs(ff_flux_[f]) / m.GetVolume(c));
      }
    }
    m.Reduce(&maxvel_, "max");
  }
  if (sem("cfl")) {
    Scal dt = var.Double["cfl"] / maxvel_;
    dt = std::min(dt, var.Double["dtmax"]);
    as_->SetTimeStep(dt);
    var.Double.Set("dt", as_->GetTimeStep());
  }
  if (sem.Nested("start")) {
    as_->StartStep();
  }
  if (sem.Nested("iter")) {
    as_->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as_->FinishStep();
  }
  if (sem("stat")) {
    sumu_ = 0.;
    auto& u = as_->GetField();
    for (auto c : m.Cells()) {
      sumu_ += u[c];
    }
    m.Reduce(&sumu_, "sum");
  }
  if (sem("t")) {
    if (IsLead()) {
      ++var.Int["iter"];
      var.Double["t"] = as_->GetTime();
    }
    if (IsRoot()) {
      auto t = var.Double["t"];
      auto dt = var.Double["dt"];
      if (!var.Double("laststat")) {
        var.Double.Set("laststat", -1e10);
      }
      if (t - var.Double["laststat"] >= var.Double["statdt"]) {
        std::cout 
            << "t=" << t 
            << " dt=" << dt 
            << std::setprecision(16) << " sumu=" << sumu_
            << std::endl;
        var.Double["laststat"] = t;
      }
    }
  }
  if (sem("comm")) { // comm for GetGlobalField
    auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
    m.Comm(&u);
  }
  sem.LoopEnd();
}

template <class M_>
class DistrSolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  using K = Advection<M>;
  using KF = KernelMeshParFactory<M, K>;
  using D = DistrMesh<KernelMeshFactory<M>>; 
  using Par = typename K::Par;

  DistrSolver(MPI_Comm comm, Vars& var, Par& par)
      : var(var), par_(par)
  {
    // Create kernel factory
    KF kf(par_);

    var.Double.Set("t", 0.);

    // Initialize blocks
    Distr* b; // base
    if (var.Int["loc"]) {
      b = CreateLocal(comm, kf, var);
    } else {
      //b = CreateCubism(comm, kf, var);
    }

    d_.reset(dynamic_cast<D*>(b));
    if (!d_) {
      throw std::runtime_error("DistrSolver: Can't cast to D");
    }
  }
  void MakeStep() {
    d_->Run();
  }
  GBlock<IdxCell, dim> GetBlock() const {
    return d_->GetGlobalBlock();
  }
  FieldCell<Scal> GetField() const {
    return d_->GetGlobalField(0);
  }
  double GetTime() const { 
    return var.Double["t"]; 
  }
  double GetTimeStep() const { 
    return var.Double["dt"]; 
  }

 private:
  Vars& var;
  std::unique_ptr<D> d_;
  Par& par_;
};

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  typename Advection<M>::Par par;
  par.fu0 = CreateInitU<Vect>(var);
  par.fvel = CreateInitVel<Vect>(var);

  DistrSolver<M> ds(comm, var, par);

  size_t k = 0;

  // Comm initial field (needed for GetField())
  ds.MakeStep();
  Dump(ds.GetField(), ds.GetBlock(), "u" + std::to_string(k) + ".dat");
  ++k;

  auto& ttmax = var.Double["ttmax"];
  auto& tmax = var.Double["tmax"];
  // Time steps until reaching ttmax
  while (tmax < ttmax) {
    std::cout 
        << "tmax = " << tmax 
        << " dt = " << ds.GetTimeStep()
        << std::endl;

    // Time steps until reaching tmax
    tmax += var.Double["dumpdt"];
    ds.MakeStep();

    Dump(ds.GetField(), ds.GetBlock(), "u" + std::to_string(k) + ".dat");
    ++k;
  }
}

int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool isroot = (!rank);

  Vars var;   // parameter storage
  Interp ip(var); // interpretor (parser)

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

  bool loc = var.Int["loc"];

  MPI_Comm comm;
  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      Main(comm, var);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, var);
  }

  MPI_Finalize();	
  return 0;
}
