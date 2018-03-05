#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>

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
#include "hydro/advection.hpp"

// Test design rules:

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

  Advection(Vars& par, const MyBlockInfo& bi);
  void Run() override;

 protected:
  using KM::par;
  using KM::m;

 private:
  void TestSolve(std::function<Scal(Vect)> fi /*initial*/,
                 std::function<Scal(Vect)> fe /*exact*/,
                 size_t cg /*check gap (separate from boundary)*/,
                 std::string name);

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;

  geom::FieldCell<Scal> fc_exact_;
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldCell<Scal> fc_src_;
  using AS = solver::AdvectionSolverExplicit<M, geom::FieldFace<Scal>>;
  std::unique_ptr<AS> as_; // advection solver
  MIdx gs_; // global mesh size
  Vect ge_; // global extent
  bool broot_; // block root
};

template <class _M>
class AdvectionFactory : public KernelMeshFactory<_M> {
 public:
  using M = _M;
  using K = Advection<M>;
  K* Make(Vars& par, const MyBlockInfo& bi) override {
    return new K(par, bi);
  }
};

template <class M>
Advection<M>::Advection(Vars& par, const MyBlockInfo& bi) 
  : KernelMesh<M>(par, bi)
{
  broot_ = (m.GetInBlockCells().GetBegin() == MIdx(0));
}

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-10;
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
#define PFCMP(a, b) \
  if (!Cmp(a, b)) { \
    std::cerr \
      << std::scientific << std::setprecision(16) \
      << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
    assert(false);\
  }


template <class M>
void Advection<M>::TestSolve(
    std::function<Scal(Vect)> fi /*initial*/,
    std::function<Scal(Vect)> fe /*exact*/,
    size_t cg /*check gap (separate from boundary)*/,
    std::string name) {
  auto sem = m.GetSem("TestSolve");
  auto& bc = m.GetBlockCells();
  if (sem("init")) {
    if (broot_) {
      std::cerr << name << std::endl;
    }
    // initial field for advection
    FieldCell<Scal> fc_u(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_u[i] = fi(x);
    }

    // zero-derivative boundary conditions for advection
    geom::MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond;
    for (auto idxface : m.Faces()) {
      if (!m.IsInner(idxface)) {
        mf_cond[idxface] =
            std::make_shared
            <solver::ConditionFaceDerivativeFixed<Scal>>(Scal(0));
      }
    }

    // velocity and flux
    const Vect vel(par.Vect["vel"]);
    ff_flux_.Reinit(m);
    for (auto idxface : m.Faces()) {
      ff_flux_[idxface] = vel.dot(m.GetSurface(idxface));
    }
    
    // source
    fc_src_.Reinit(m, 0.);

    as_.reset(new AS(m, fc_u, mf_cond, &ff_flux_, 
              &fc_src_, 0., par.Double["dt"]));

    // exact solution
    fc_exact_.Reinit(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_exact_[i] = fe(x);
    }
  }
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
    PFCMP(Mean(fc_exact_, m, mask), Mean(fc, m, mask));
    PFCMP(DiffMax(fc_exact_, fc, m, mask), 0.);
  }
}

template <class M>
void Advection<M>::Run() {
  par.Double["extent"] = 1.; // TODO don't overwrite extent
  Scal extent = par.Double["extent"];
  {
    MIdx p(par.Int["px"], par.Int["py"], par.Int["pz"]);
    MIdx b(par.Int["bx"], par.Int["by"], par.Int["bz"]);
    using B = MyBlock;
    MIdx s(B::sx, B::sy, B::sz); // block size inner
    gs_ = p * b * s;
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
  auto f0 = [](Vect x) { return 0.; };
  auto f1 = [](Vect x) { return 1.; };
  auto flin = [](Vect x) { 
      Vect b(0.5); 
      x -= b;
      Scal a = x[0] + 0.3 * x[1] - 0.5 * x[2]; 
      return a > 0. ? 1. : 0.;
      //return std::sin(a);
    };
  auto fline = [=](Vect x) { 
      x -= vel * t; 
      return flin(x);
  };

  if (sem.Nested()) {
    TestSolve(f0, f0, 0, "f=0");
  }
  if (sem.Nested()) {
    TestSolve(f1, f1, 0, "f=1");
  }
  if (sem.Nested()) {
    TestSolve(flin, fline, nc, "f=step");
  }
}


void Main(MPI_Comm comm, bool loc, Vars& par) {
  // read config files, parse arguments, maybe init global fields
  using M = geom::MeshStructured<Scal, 3>;
  using K = Advection<M>;
  using KF = AdvectionFactory<M>;
  using D = Distr;
  
  KF kf;

  const int es = 8;
  const int hl = par.Int["hl"];
  const int bs = 16;
  
  // Initialize buffer mesh and make kernel for each block.
  std::unique_ptr<Distr> d;
  if (loc) {
    d = CreateLocal(comm, kf, bs, es, hl, par);
  } else {
    d = CreateCubism(comm, kf, bs, es, hl, par);
  }

  d->Step();
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
