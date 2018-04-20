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
#include "hydro/advection_vof.hpp"
#include "hydro/advection_tvd.hpp"

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
            std::function<Vect(Vect,Scal)> fvel /*vel(x,t)*/);
  void Run() override;

 protected:
  using KM::par;
  using KM::m;
  using KM::IsRoot;
  using KM::IsLead;

 private:
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldCell<Scal> fc_src_;
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

template <class _M>
class AdvectionFactory : public KernelMeshFactory<_M> {
 public:
  using M = _M;
  using K = Advection<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  AdvectionFactory(std::function<Scal(Vect)> fu0,
                   std::function<Vect(Vect,Scal)> fvel) 
      : fu0_(fu0), fvel_(fvel) {}
  K* Make(Vars& par, const MyBlockInfo& bi) override {
    return new K(par, bi, fu0_, fvel_);
  }

 private:
  std::function<Scal(Vect)> fu0_;
  std::function<Vect(Vect,Scal)> fvel_;
};

template <class M>
Advection<M>::Advection(Vars& par, const MyBlockInfo& bi,
                        std::function<Scal(Vect)> fu0,
                        std::function<Vect(Vect,Scal)> fvel /*(x,t)*/)
    : KernelMesh<M>(par, bi)
    , fvel_(fvel)
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
    ff_flux_[f] = fvel_(x, 0.).dot(m.GetSurface(f));
  }
  
  // source
  fc_src_.Reinit(m, 0.);

  std::string as = par.String["advection_solver"];
  if (as == "tvd") {
    auto& p = aspar1_;
    p.sharp = par.Double["sharp"];
    p.sharpo = par.Double["sharpo"];
    p.split = par.Int["split"];
    as_.reset(new AS1(m, u0, bc, &ff_flux_, 
              &fc_src_, 0., par.Double["dt"], &p));
  } else if (as == "vof") {
    auto& p = aspar2_;
    as_.reset(new AS2(m, u0, bc, &ff_flux_, 
              &fc_src_, 0., par.Double["dt"], &p));
  } else {
    throw std::runtime_error("Unknown advection_solver=" + as);
  }
}


template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("advection");
  if (sem("comm")) { // comm for GetGlobalField, initial
    auto& u = const_cast<geom::FieldCell<Scal>&>(as_->GetField());
    m.Comm(&u);
  }
  sem.LoopBegin();
  if (sem("empty")) {
    //nop // TODO: bugfix loop, empty stage needed
  }
  if (sem("checkloop")) {
    if (as_->GetTime() >= par.Double["tmax"]) {
      sem.LoopBreak();
    }
  }
  if (sem("vel")) {
    const Scal t = as_->GetTime();
    for (auto f : m.AllFaces()) {
      Vect x = m.GetCenter(f);
      ff_flux_[f] = fvel_(x, t).dot(m.GetSurface(f));
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
    Scal dt = par.Double["cfl"] / maxvel_;
    dt = std::min(dt, par.Double["dtmax"]);
    as_->SetTimeStep(dt);
    par.Double.Set("dt", as_->GetTimeStep());
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
      ++par.Int["iter"];
      par.Double["t"] = as_->GetTime();
    }
    if (IsRoot()) {
      auto t = par.Double["t"];
      auto dt = par.Double["dt"];
      if (!par.Double("laststat")) {
        par.Double.Set("laststat", -1e10);
      }
      if (t - par.Double["laststat"] >= par.Double["statdt"]) {
        std::cout 
            << "t=" << t 
            << " dt=" << dt 
            << std::setprecision(16) << " sumu=" << sumu_
            << std::endl;
        par.Double["laststat"] = t;
      }
    }
  }
  if (sem("comm")) { // comm for GetGlobalField
    auto& u = const_cast<geom::FieldCell<Scal>&>(as_->GetField());
    m.Comm(&u);
  }
  sem.LoopEnd();
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
void Dump(const geom::FieldCell<Scal>& u, const B& b, std::string on) {
  std::ofstream o(on.c_str());

  auto s = b.GetDimensions();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : b.Range()) {
    o << u[c] << " ";
  }

  o << std::endl;
}

template <class M_>
class DistrSolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using IdxCell = geom::IdxCell;

  using KF = AdvectionFactory<M>;
  using D = DistrMesh<KernelMeshFactory<M>>; 

  DistrSolver(MPI_Comm comm, Vars& apar, 
              std::function<Scal(Vect)> fu0,
              std::function<Vect(Vect,Scal)> fvel /*vel(x,t)*/)
      : par(apar)
  {
    // Create kernel factory
    KF kf(fu0, fvel);

    par.Double.Set("t", 0.);

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
  void MakeStep() {
    d_->Run();
  }
  geom::GBlock<IdxCell, dim> GetBlock() const {
    return d_->GetGlobalBlock();
  }
  geom::FieldCell<Scal> GetField() const {
    return d_->GetGlobalField(0);
  }
  double GetTime() const { 
    return par.Double["t"]; 
  }
  double GetTimeStep() const { 
    return par.Double["dt"]; 
  }

 private:
  Vars& par;
  std::unique_ptr<D> d_;
};

void Main(MPI_Comm comm, Vars& par) {
  using M = geom::MeshStructured<double, 3>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  std::function<Scal(Vect)> fu0;
  {
    std::string v = par.String["init_vf"];
    if (v == "circle") {
      auto fucircle = [](Vect x) -> Scal { 
        return Vect(0.5, 0.263662, 0.).dist(x) < 0.2 ? 1. : 0.; 
      };
      fu0 = fucircle;
    } else if (v == "sinc") {
      Vect k(par.Vect["sinck"]);
      auto fusinc = [k](Vect x) -> Scal { 
        x -= Vect(0.5);
        x *= k;
        Scal r = x.norm();
        Scal u0 = -0.2;
        Scal u1 = 1.;
        Scal u = std::sin(r) / r;
        u = (u - u0) / (u1 - u0);
        return std::max(0., std::min(1., u));
      };
      fu0 = fusinc;
    } else {
      throw std::runtime_error("Unknown init_vf=" + v);
    }
  }


  std::function<Vect(Vect,Scal)> fvel;
  {
    std::string v = par.String["init_vel"];
    if (v == "uni") {
      Vect vel(par.Vect["vel"]);
      auto fveluni = [vel](Vect x, Scal /*t*/) -> Vect { 
        return Vect(vel); 
      };
      fvel = fveluni;
    } else if (v == "sincos") {
      Scal revt = par.Double["revt"]; // reverse time
      auto fvelsincos = [revt](Vect x, Scal t) -> Vect { 
        x = x * M_PI;
        Vect r(std::sin(x[0]) * std::cos(x[1]), 
               -std::cos(x[0]) * std::sin(x[1]),
               0.);
        if (t > revt) { r *= -1.; }
        return r; 
      };
      fvel = fvelsincos;
    } else {
      throw std::runtime_error("Unknown init_vel=" + v);
    }
  }

  DistrSolver<M> ds(comm, par, fu0, fvel);

  size_t k = 0;

  // Comm initial field (needed for GetField())
  ds.MakeStep();
  Dump(ds.GetField(), ds.GetBlock(), "u" + std::to_string(k) + ".dat");
  ++k;

  auto& ttmax = par.Double["ttmax"];
  auto& tmax = par.Double["tmax"];
  // Time steps until reaching ttmax
  while (tmax < ttmax) {
    std::cout 
        << "tmax = " << tmax 
        << " dt = " << ds.GetTimeStep()
        << std::endl;

    // Time steps until reaching tmax
    tmax += par.Double["dumpdt"];
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
