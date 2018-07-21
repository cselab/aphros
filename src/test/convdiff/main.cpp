#undef NDEBUG
#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>
#include <utility>
#include <tuple>
#include <sstream>
#include <memory>

#include "parse/vars.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
#include "util/suspender.h"
#include "geom/vect.h"
#include "geom/mesh.h"
#include "linear/linear.h"
#include "solver/solver.h"
#include "solver/convdiffi.h"
#include "dump/output.h"
#include "dump/output_paraview.h"

struct GPar {};

template <class M_>
class Convdiff : public KernelMeshPar<M_, GPar> {
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
  using P::m;
  using P::IsRoot;

 private:
  void TestSolve(std::function<Scal(Vect)> fi /*initial*/,
                 std::function<Scal(Vect)> fe /*exact*/,
                 size_t cg /*check gap (separate from boundary)*/,
                 std::string name,
                 bool check /*abort if differs from exact*/);

  template <class T>
  using FieldCell = FieldCell<T>;
  template <class T>
  using FieldFace = FieldFace<T>;

  FieldCell<Scal> fc_exact_;
  FieldFace<Scal> ff_flux_;
  FieldCell<Scal> fc_src_;
  using AS = solver::ConvectionDiffusionScalarImplicit<M>;
  std::unique_ptr<AS> as_;
  MIdx gs_; // global mesh size
  Vect ge_; // global extent
  FieldCell<Scal> fc_; // buffer
  FieldCell<Scal> fc_sc_; // scaling
  FieldFace<Scal> ff_d_; // diffusion rate
};

template <class T>
bool Cmp(T a, T b) {
  return a == b;
}

template <>
bool Cmp<double>(double a, double b) {
  return std::abs(a - b) < 1e-6;
}


template <class Idx, class B, class Scal>
Scal DiffMax(
    const B& b,
    const GField<Scal, Idx>& u,
    const GField<Scal, Idx>& v,
    const GField<bool, Idx>& mask) {
  Scal r = 0;
  for (auto i : GRange<Idx>(b)) {
    if (mask[i]) {
      r = std::max(r, std::abs(u[i] - v[i]));
    }
  }
  return r;
}

template <class Idx, class B, class Scal>
Scal Max(
    const B& b,
    const GField<Scal, Idx>& u,
    const GField<bool, Idx>& mask) {
  Scal r = 0;
  for (auto i : GRange<Idx>(b)) {
    if (mask[i]) {
      r = std::max(r, u[i]);
    }
  }
  return r;
}

template <class Idx, class B, class Scal>
Scal Mean(
    const B& b,
    const GField<Scal, Idx>& u,
    const GField<bool, Idx>& mask) {
  Scal r = 0;
  Scal w = 0.;
  for (auto i : GRange<Idx>(b)) {
    if (mask[i]) {
      r += u[i];
      w += 1.;
    }
  }
  return r / w;
}


template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v,
    const M& m,
    const GField<bool, Idx>& mask) {
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
    const GField<typename M::Scal, Idx>& u,
    const M& m,
    const GField<bool, Idx>& mask) {
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
    const GField<typename M::Scal, Idx>& u,
    const M& m,
    const GField<bool, Idx>& mask) {
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
#define PCMP(a, b, fatal) \
  std::cerr \
    << std::scientific << std::setprecision(16) \
    << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  if (!Cmp(a, b)) { \
    assert(!fatal);\
  }


// Print CMP if false
#define PFCMP(a, b, fatal) \
  if (!Cmp(a, b)) { \
    std::cerr \
      << std::scientific << std::setprecision(16) \
      << "Failed cmp: " << std::endl \
      << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
    assert(!fatal);\
  }


template <class M>
void Convdiff<M>::TestSolve(
    std::function<Scal(Vect)> fi /*initial*/,
    std::function<Scal(Vect)> fe /*exact*/,
    size_t cg /*check gap (separate from boundary)*/,
    std::string name,
    bool check /*abort if different from exact*/) {
  auto sem = m.GetSem("TestSolve");
  auto& bc = m.GetBlockCells();
  if (sem("init")) {
    if (IsRoot()) {
      std::cerr << name << std::endl;
    }
    // initial field 
    FieldCell<Scal> fc_u(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_u[i] = fi(x);
    }

    // global mesh size
    MIdx gs;
    {
      MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
      MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
      MIdx bs(var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]);
      gs = p * b * bs;
    }

    using Dir = typename M::Dir;

    MapFace<std::shared_ptr<solver::CondFace>> mf_cond;

    // boundary xm of global mesh
    auto gxm = [this](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::i &&
          m.GetBlockFaces().GetMIdx(i)[0] == 0;
    };
    auto gxp = [this,gs](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::i &&
          m.GetBlockFaces().GetMIdx(i)[0] == gs[0];
    };
    auto gym = [this](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::j &&
          m.GetBlockFaces().GetMIdx(i)[1] == 0;
    };
    auto gyp = [this,gs](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::j &&
          m.GetBlockFaces().GetMIdx(i)[1] == gs[1];
    };
    auto gzm = [this](IdxFace i) -> bool {
      return dim >= 3 && m.GetDir(i) == Dir::k &&
          m.GetBlockFaces().GetMIdx(i)[2] == 0;
    };
    auto gzp = [this,gs](IdxFace i) -> bool {
      return dim >= 3 && m.GetDir(i) == Dir::k &&
          m.GetBlockFaces().GetMIdx(i)[2] == gs[2];
    };
    auto parse = [](std::string s, IdxFace, size_t nci, M&) 
       -> std::shared_ptr<solver::CondFace> {
      std::stringstream arg(s);

      std::string name;
      arg >> name;

      if (name == "value") {
        Scal a;
        arg >> a;
        return std::make_shared <solver::
            CondFaceValFixed<Scal>>(a, nci);
      } else if (name == "derivative") {
        Scal a;
        arg >> a;
        return std::make_shared <solver::
            CondFaceGradFixed<Scal>>(a, nci);
      } else {
        assert(false);
      }
    };
    // Set condition bc for face i on global box boundary
    // choosing proper neighbour cell id (nci)
    // Return true if on global boundary
    auto set_bc = [&](IdxFace i, std::string bc) -> bool {
      if (gxm(i) || gym(i) || gzm(i)) {
        mf_cond[i] = parse(bc, i, 1, m);
        return true;
      } else if (gxp(i) || gyp(i) || gzp(i)) {
        mf_cond[i] = parse(bc, i, 0, m);
        return true;
      } 
      return false;
    };

    // Boundary conditions for fluid 
    auto ff = m.AllFaces();
    std::vector<std::pair<std::string, std::function<bool(IdxFace)>>> pp = 
        {{"bc_xm", gxm}, {"bc_xp", gxp},
         {"bc_ym", gym}, {"bc_yp", gyp},
         {"bc_zm", gzm}, {"bc_zp", gzp}};

    for (auto p : pp) {
      if (auto bc = var.String(p.first)) {
        for (auto i : ff) {
          p.second(i) && set_bc(i, *bc);
        }
      } 
    }

    // velocity and flux
    const Vect vel(var.Vect["vel"]);
    ff_flux_.Reinit(m);
    for (auto idxface : m.AllFaces()) {
      ff_flux_[idxface] = vel.dot(m.GetSurface(idxface));
    }

    // cell conditions 
    // (empty)
    MapCell<std::shared_ptr<solver::CondCell>> mc_cond;
    
    // source
    fc_src_.Reinit(m, 0.);

    // diffusion rate
    ff_d_.Reinit(m, var.Double["mu"]);

    // scaling for as_
    fc_sc_.Reinit(m, 1.); 

    auto p = std::make_shared<typename AS::Par>();
    p->relax = var.Double["relax"];
    p->guessextra = 0.;
    p->second = 0;

    as_.reset(new AS(
          m, fc_u, mf_cond, mc_cond, 
          &fc_sc_, &ff_d_, &fc_src_, &ff_flux_,
          0., var.Double["dt"], p));

    // exact solution
    fc_exact_.Reinit(m);
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_exact_[i] = fe(x);
    }
  }
  for (size_t n = 0; n < size_t(var.Int["num_steps"]); ++n) {
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
    GBlockCells<dim> cbc(MIdx(cg), gs_ - MIdx(2 * cg)); // check block
    FieldCell<bool> mask(m, false);
    for (auto i : m.AllCells()) {
      if (cbc.IsInside(bc.GetMIdx(i))) {
        mask[i] = true;
      }
    }
    auto& fc = as_->GetField();
    if (check) {
      PCMP(Mean(fc_exact_, m, mask), Mean(fc, m, mask), check);
      PCMP(DiffMax(fc_exact_, fc, m, mask), 0., check);
    }
  }
  if (sem("comm")) {
    fc_ = as_->GetField();
    m.Comm(&fc_);
    m.Comm(&fc_exact_);
  }
}

template <class M>
void Convdiff<M>::Run() {
  var.Double["extent"] = 1.; // TODO don't overwrite extent
  Scal extent = var.Double["extent"];
  {
    MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
    MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
    MIdx bs(var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]);
    gs_ = p * b * bs;
  }
  Scal dx = extent / gs_.norminf(); 
  assert(dx > 0. && dx < extent);
  ge_ = Vect(gs_) * dx;
  Vect vel(var.Vect["vel"]);
  Scal ns = var.Int["num_steps"];
  if (auto pcfl = var.Double("cfl")) {
    Scal cfl = *pcfl;
    Scal dt = dx * cfl / vel.norminf();
    var.Double.Set("dt", dt);
    //Scal nc = cfl * ns; // distance in cells passed
  }
  Scal t = var.Double["dt"] * ns;

  auto sem = m.GetSem("Run");
  auto f = [](Vect x) { 
      return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]); 
    };
  //auto f0 = [](Vect { return 0.; };
  //auto f1 = [](Vect) { return 1.; };
  auto fx = [](Vect x) { return Vect(1., -1., 1.).dot(x); };
  auto fex = [=](Vect x) { x -= vel * t; return fx(x); };

  if (sem.Nested()) {
    TestSolve(fx, fex, 6, "fx", false);
  }
  if (sem.Nested()) {
    TestSolve(f, f, 0, "f", false);
  }
}

using Scal = double;
using M = MeshStructured<Scal, 3>;
using K = Convdiff<M>;
using BC = typename M::BlockCells;
using FC = FieldCell<Scal>;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

std::tuple<BC, FC, FC> Solve(MPI_Comm comm, Vars& var) {
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();

  return std::make_tuple(ds.GetBlock(), ds.GetField(0), ds.GetField(1));
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

void Main(MPI_Comm comm, Vars& var0) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  bool isroot = (!rank);

  Vars var = var0;
  
  if (isroot) {
    std::cerr << "solve ref" << std::endl;
  }

  BC b, bb;
  FC fa, fea, fb, feb;
  std::tie(b, fa, fea) = Solve(comm, var); // fa non-empty only on root

  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx ms = b.GetDimensions();
  auto m = InitUniformMesh<M>(dom, MIdx(0), ms, 0, true, ms);

  if (isroot) {
    std::cerr << "solve bs/2" << std::endl;
  }

  var.Int["bsx"] /= 2;
  var.Int["bsy"] /= 2;
  var.Int["bsz"] /= 2;
  var.Int["bx"] *= 2;
  var.Int["by"] *= 2;
  var.Int["bz"] *= 2;
  std::tie(bb, fb, feb) = Solve(comm, var); // fb non-empty only on root

  if (isroot) {
    PCMP(b.GetEnd(), bb.GetEnd(), true);

    FieldCell<bool> mask(b, true);

    Dump({&fa, &fea, &fb}, {"a", "ea", "b"}, m, "a");

    PCMP(Mean(b, fa, mask), Mean(b, fb, mask), true);
    PCMP(DiffMax(b, fa, fb, mask), 0., true);
  }
}


int main (int argc, const char ** argv) {
  return RunMpi(argc, argv, Main);
}
