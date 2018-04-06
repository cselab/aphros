include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>
#include <utility>
#include <tuple>

#include "CubismDistr/Vars.h"
#include "CubismDistr/Interp.h"
#include "CubismDistr/Kernel.h"
#include "CubismDistr/KernelMesh.h"
#include "CubismDistr/ICubism.h"
#include "CubismDistr/ILocal.h"
#include "CubismDistr/Distr.h"

#include "hydro/suspender.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/linear.hpp"
#include "hydro/advection.hpp"

#include "hydro/output.hpp"
#include "hydro/output_paraview.hpp"


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
  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;

  geom::FieldCell<Scal> fc_exact_;
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldCell<Scal> fc_src_;
  using AS = solver::AdvectionSolverExplicit<M>;
  std::unique_ptr<AS> as_; // advection solver
  MIdx gs_; // global mesh size
  Vect ge_; // global extent
  geom::FieldCell<Scal> fc_; // buffer
  FieldCell<Scal> fc_sc_; // scaling
  FieldFace<Scal> ff_d_; // diffusion rate
};

template <class _M>
class AdvectionFactory : public KernelMeshFactory<_M> {
 public:
  using M = _M;
  using K = Advection<M>;
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
  FieldCell<Scal> u0(m, 0);
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

  // diffusion rate
  ff_d_.Reinit(m, par.Double["mu"]);

  as_.reset(new AS(m, u0, bc, &ff_flux_, 
            &fc_src_, 0., par.Double["dt"], 
            0., 0., 1.));
}

template <class T>
bool Cmp(T a, T b) {
  return a == b;
}

template <>
bool Cmp<double>(double a, double b) {
  return std::abs(a - b) < 1e-10;
}


template <class Idx, class B, class Scal>
Scal DiffMax(
    const B& b,
    const geom::GField<Scal, Idx>& u,
    const geom::GField<Scal, Idx>& v,
    const geom::GField<bool, Idx>& mask) {
  Scal r = 0;
  for (auto i : geom::GRange<Idx>(b)) {
    if (mask[i]) {
      r = std::max(r, std::abs(u[i] - v[i]));
    }
  }
  return r;
}

template <class Idx, class B, class Scal>
Scal Max(
    const B& b,
    const geom::GField<Scal, Idx>& u,
    const geom::GField<bool, Idx>& mask) {
  Scal r = 0;
  for (auto i : geom::GRange<Idx>(b)) {
    if (mask[i]) {
      r = std::max(r, u[i]);
    }
  }
  return r;
}

template <class Idx, class B, class Scal>
Scal Mean(
    const B& b,
    const geom::GField<Scal, Idx>& u,
    const geom::GField<bool, Idx>& mask) {
  Scal r = 0;
  Scal w = 0.;
  for (auto i : geom::GRange<Idx>(b)) {
    if (mask[i]) {
      r += u[i];
      w += 1.;
    }
  }
  return r / w;
}


template <class Idx, class M>
typename M::Scal DiffMax(
    const geom::GField<typename M::Scal, Idx>& u,
    const geom::GField<typename M::Scal, Idx>& v,
    const M& m,
    const geom::GField<bool, Idx>& mask) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    if (mask[i]) {
      r = std::max(r, std::abs(u[i] - v[i]));
    }
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Max(
    const geom::GField<typename M::Scal, Idx>& u,
    const M& m,
    const geom::GField<bool, Idx>& mask) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    if (mask[i]) {
      r = std::max(r, u[i]);
    }
  }
  return r;
}


template <class Idx, class M>
typename M::Scal Mean(
    const geom::GField<typename M::Scal, Idx>& u,
    const M& m,
    const geom::GField<bool, Idx>& mask) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0.;
  for (auto i : m.template Get<Idx>()) {
    if (mask[i]) {
      r += u[i];
      w += 1.;
    }
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


// Print CMP if false
#define PFCMP(a, b, ftl) \
  if (!Cmp(a, b)) { \
    std::cerr \
      << std::scientific << std::setprecision(16) \
      << "Failed cmp: " << std::endl \
      << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
    assert(!ftl);\
  }


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

template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("TestSolve");
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
using KF = AdvectionFactory<M>;
using D = DistrMesh<KernelMeshFactory<M>>;
using BC = typename M::BlockCells;
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

// Solver derived from UnsteadySolver
// Constructor:
//   initial field Scal(Vect x)
//   velocity Vect(Vect x)
//   minimal par (bx, bpx, px, extent, dt, tmax)
//   no config files
// GBlock<IdxCell, 3> GetBlock() // global block
// FieldCell<Scal> GetField()    // current field on global block
//
// All necessary communication is done internally.
// Maybe restrict some calls to root.

// Run distributed solver and return result on global mesh
std::tuple<BC, FC, FC> Solve(MPI_Comm comm, Vars& par) {
  KF kf;

  bool loc = par.Int["loc"];
  
  Distr* dr;
  if (loc) {
    dr = CreateLocal(comm, kf, par);
  } else {
    dr = CreateCubism(comm, kf, par);
  }

  std::unique_ptr<D> d(dynamic_cast<D*>(dr));
  assert(d);
  d->Run();

  return std::make_tuple(
      d->GetGlobalBlock(), d->GetGlobalField(0), d->GetGlobalField(1));
}

void Dump(std::vector<const FC*> u, std::vector<std::string> un, 
          M& m, std::string fn="a") {
  output::Content con;
  for (size_t k = 0; k < u.size(); ++k) {
    auto p = u[k];
    con.emplace_back(
        new output::EntryFunction<Scal, IdxCell, M>(
            un[k], m, [p](IdxCell i) { return (*p)[i]; }));
  }

  output::SessionParaviewStructured<M> ses(con, "", fn /*filename*/, m);
  ses.Write(0, "t");
}

class DistrSolver : public solver::UnsteadySolver {
 public:
  using Scal = double;
  using AS = solver::AdvectionSolverExplicit<M>;
  using Mesh = geom::MeshStructured<Scal, 3>;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using IdxCell = geom::IdxCell;
  DistrSolver(MPI_Comm comm, Vars par, 
              std::function<Scal(Vect)> fu0,
              std::function<Vect(Vect)> fvel)
      : comm_(comm), par(par)
  {

  }

 private:
  MPI_Comm comm_;
  Vars par;
};

void Main(MPI_Comm comm, Vars& par0) {
  Vars par = par0;
  BC b;
  FC f, fe;
  par.Int["bsx"] = 32;
  par.Int["bsy"] = 32;
  par.Int["bsz"] = 1;
  par.Int["bx"] = 1;
  par.Int["by"] = 1;
  par.Int["bz"] = 1;
  par.Int["px"] = 1;
  par.Int["py"] = 1;
  par.Int["pz"] = 1;
  par.Int["fatal"] = 0;
  Dump(f, b, "u0.dat");
  std::tie(b, f, fe) = Solve(comm, par);
  Dump(f, b, "u.dat");

  return;
  
  std::cerr << "solve with bsx=" << par.Int["bsx"] << std::endl;
  BC b0, b1;
  FC f0, fe0, f1, fe1;
  std::tie(b0, f0, fe0) = Solve(comm, par);


  // Create uniform mesh in unit cube
  geom::Rect<Vect> dom(Vect(0), Vect(1));
  auto m = geom::InitUniformMesh<M>(dom, MIdx(0), b0.GetDimensions(), 0);

  par.Int["bsx"] /= 2;
  par.Int["bsy"] /= 2;
  par.Int["bsz"] /= 2;
  par.Int["bx"] *= 2;
  par.Int["by"] *= 2;
  par.Int["bz"] *= 2;
  std::cerr << "solve with bsx=" << par.Int["bsx"] << std::endl;
  std::tie(b1, f1, fe1) = Solve(comm, par);

  // Check both grids have same size (assume Begin=0)
  PCMP(b0.GetEnd(), b1.GetEnd());

  geom::FieldCell<bool> mask(b0, true);

  Dump({&f0, &fe0, &f1}, {"0", "0e", "1"}, m, "a");

  PCMP(Mean(b0, f0, mask), Mean(b1, f1, mask));
  PCMP(DiffMax(b0, f0, f1, mask), 0.);
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
