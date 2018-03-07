#include <iostream>
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
                 std::string name,
                 bool check /*abort if differs from exact*/);

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
  geom::FieldCell<Scal> fc_; // buffer
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

template <class T>
bool Cmp(T a, T b) {
  return a == b;
}

template <>
bool Cmp<Scal>(Scal a, Scal b) {
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
    std::string name,
    bool check /*abort if different from exact*/) {
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
    if(0)
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
    for (auto idxface : m.AllFaces()) {
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
  if (check && sem("check")) {
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
  if (sem("comm")) {
    fc_ = as_->GetField();
    m.Comm(&fc_);
  }
}

template <class M>
void Advection<M>::Run() {
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
  auto fx = [](Vect x) { return x[0]; };
  auto fex = [=](Vect x) { x -= vel * t; return fx(x); };

  if (sem.Nested()) {
    TestSolve(f0, f0, 0, "f0", true);
  }
  if (sem.Nested()) {
    TestSolve(f1, f1, 0, "f1", true);
  }
  if (sem.Nested()) {
    TestSolve(fx, fex, 3, "fx", false);
  }
  if (sem.Nested()) {
    //TestSolve(f, f, 0, "f", false);
  }
}

using M = geom::MeshStructured<Scal, 3>;
using K = Advection<M>;
using KF = AdvectionFactory<M>;
using D = DistrMesh<KernelMeshFactory<M>>;
using BC = typename M::BlockCells;
using FC = geom::FieldCell<Scal>;
using IdxCell = geom::IdxCell;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

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
  d->Step();

  return std::make_tuple(
      d->GetGlobalBlock(), d->GetGlobalField(0), d->GetGlobalField(0));
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

void Main(MPI_Comm comm, Vars& par0) {
  Vars par = par0;
  
  std::cerr << "solve ref" << std::endl;
  BC b, bb;
  FC fa, fea, fb, feb;
  std::tie(b, fa, fea) = Solve(comm, par);

  using Scal = double;
  geom::Rect<Vect> dom(Vect(0), Vect(1));
  auto m = geom::InitUniformMesh<M>(dom, MIdx(0), b.GetDimensions(), 0);

  std::cerr << "solve bs/2" << std::endl;
  par.Int["bsx"] /= 2;
  par.Int["bsy"] /= 2;
  par.Int["bsz"] /= 2;
  par.Int["bx"] *= 2;
  par.Int["by"] *= 2;
  par.Int["bz"] *= 2;
  std::tie(bb, fb, feb) = Solve(comm, par);

  PCMP(b.GetEnd(), bb.GetEnd());

  geom::FieldCell<bool> mask(b, true);

  Dump({&fa, &fea}, {"u", "ue"}, m);

  PCMP(Mean(b, fa, mask), Mean(b, fb, mask));
  PCMP(DiffMax(b, fa, fb, mask), 0.);
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
      Main(comm, par);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, par);
  }


  MPI_Finalize();	
  return 0;
}
