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

#include "parse/vars.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"

#include "util/suspender.h"
#include "geom/vect.h"
#include "geom/mesh.h"
#include "solver/solver.h"
#include "solver/advection.h"
#include "solver/vof.h"
#include "solver/tvd.h"
#include "func/init_vel.h"
#include "func/init_u.h"
#include "dump/dumper.h"
#include "dump/output.h"

#include "cmp.h"

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

template <class M_>
class Advection;

template <class M_>
struct GPar {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  std::function<void(FieldCell<Scal>&, const M&)> fu0; // init vf
  std::function<Vect(Vect /*x*/,Scal /*t*/)> fv; // velocity

  using K = Advection<M>;
  DistrSolver<M, K>* ds;
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
  using Sem = typename Mesh::Sem;

  //using P::P; // inherit constructor
  Advection(Vars& var, const MyBlockInfo& bi, Par& par);
  void Run() override;
  void Init(Sem& sem);
  void Dump(Sem& sem);

 protected:
  using P::var;
  using P::par_;
  using P::m;
  using P::IsRoot;
  using P::IsLead;

 private:
  FieldFace<Scal> ff_flux_;
  FieldCell<Scal> fc_src_;
  FieldCell<Scal> fcu_; // volume fraction (used for initial)
  std::unique_ptr<solver::AdvectionSolver<M>> as_; // advection solver
  Scal maxvel_; // maximum velocity relative to cell length [1/time]
                // cfl = dt * maxvel
  Scal sumu_; // sum of fluid volume
  FieldCell<Scal> fcnx_, fcny_, fcnz_; // normal to interface (tmp)
                                             // used for Vof dump
  Dumper dmf_; // fields
  Dumper dms_; // statistics
};

template <class M>
Advection<M>::Advection(Vars& v, const MyBlockInfo& b, Par& p)
    : KernelMeshPar<M, Par>(v, b, p)
    , dmf_(v, "dump_field_"), dms_(v, "dump_stat_") {}

template <class M>
void Advection<M>::Init(Sem& sem) {
  if (sem("init-local")) {
    // initial field for advection
    fcu_.Reinit(m);
    par_.fu0(fcu_, m);
    m.Comm(&fcu_);
  }

  if (sem("init-create")) {
    // flux
    ff_flux_.Reinit(m, 0);
    for (auto f : m.AllFaces()) {
      Vect x = m.GetCenter(f);
      ff_flux_[f] = par_.fv(x, 0.).dot(m.GetSurface(f));
    }

    // boundary conditions for advection (empty)
    MapFace<std::shared_ptr<solver::CondFace>> bc;

    // cell conditions for advection (empty)
    MapCell<std::shared_ptr<solver::CondCell>> mc_cond;

    // source
    fc_src_.Reinit(m, 0.);

    std::string as = var.String["advection_solver"];
    if (as == "tvd") {
      using AS = solver::AdvectionSolverExplicit<M>;
      auto p = std::make_shared<typename AS::Par>();
      p->sharp = var.Double["sharp"];
      p->sharpo = var.Double["sharpo"];
      p->split = var.Int["split"];
      as_.reset(new AS(m, fcu_, bc, &ff_flux_, 
                       &fc_src_, 0., var.Double["dt"], p));
    } else if (as == "vof") {
      using AS = solver::Vof<M>;
      auto p = std::make_shared<typename AS::Par>();
      as_.reset(new AS(m, fcu_, bc, &ff_flux_, 
                       &fc_src_, 0., var.Double["dt"], p));
    } else {
      throw std::runtime_error("Unknown advection_solver=" + as);
    }
  }
}

// Writes legacy vtk polydata 
// xx: points
// pp: polygons as lists of indices
// fn: path
// cm: comment
// Binary format.
template <class Vect>
void WriteVtkPoly(const std::vector<Vect>& xx, 
                  const std::vector<std::vector<size_t>>& pp,  
                  const std::string& fn,
                  const std::string& cm="") {
  std::ofstream f(fn.c_str());
  f << "# vtk DataFile Version 2.0\n";
  f << cm << "\n";
  f << "ASCII\n";
  f << "DATASET POLYDATA\n";

  f << "POINTS " <<  xx.size() << " float\n";
  for (auto& x : xx) {
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }

  size_t np = 0; // total number of vortices
  for (auto& p : pp) {
    np += p.size();
  }
  f << "POLYGONS " << pp.size() << " " << (np + pp.size()) << "\n";
  for (auto& p : pp) {
    f << p.size() << " " << p << "\n";
  }
}

// Converts to index representation.
// vv: polygons as lists of points
// Returns:
// xx: points
// pp: polygons as lists of indices
template <class Vect>
void Convert(const std::vector<std::vector<Vect>>& vv, 
             std::vector<Vect>& xx, 
             std::vector<std::vector<size_t>>& pp) {
  xx.resize(0);
  pp.resize(0);
  for (auto& v : vv) {
    pp.emplace_back();
    for (auto& x : v) {
      pp.back().push_back(xx.size());
      xx.push_back(x);
    }
  }
}

// Writes legacy vtk polydata 
// vv: polygons as lists of points
// fn: path
// cm: comment
// Binary format.
template <class Vect>
void WriteVtkPoly(const std::vector<std::vector<Vect>>& vv,  
                  const std::string& fn,
                  const std::string& cm="") {
  std::vector<Vect> xx;
  std::vector<std::vector<size_t>> pp;
  Convert(vv, xx, pp);
  WriteVtkPoly(xx, pp, fn, cm);
}

template <class M, class Vect=typename M::Vect>
Vect GetCellSize(const M& m) {
  Vect h; // result
  IdxCell c0(0);
  h = m.GetNode(m.GetNeighbourNode(c0, 7)) - 
      m.GetNode(m.GetNeighbourNode(c0, 0));
  assert(std::abs(h.prod() - m.GetVolume(c0)) < 1e-10);
  return h;
}

template <class M>
void Advection<M>::Dump(Sem& sem) {
  if (sem("dump")) {
    if (dmf_.Try(var.Double["t"], var.Double["dt"])) {
      auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
      m.Dump(&u, "u");
      auto& k = const_cast<FieldCell<Scal>&>(as_->GetCurv());
      m.Dump(&k, "k");
      if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
        auto& a = const_cast<FieldCell<Scal>&>(as->GetAlpha());
        m.Dump(&a, "a");
        auto &n = as->GetNormal();
        m.Dump(&n, 0, "nx");
        m.Dump(&n, 1, "ny");
        m.Dump(&n, 2, "nz");
      }

      if (IsRoot()) {
        dmf_.Report();
      }
    }
  }
  if (sem("dumpwrite")) {
    // Empty stage for DumpWrite
    // TODO: revise
  }
  // Dump reconstructed interface
  if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
    if (sem("dumpsurf")) {
      if (IsLead()) {
        auto bc = par_.ds->GetBlock();

        auto u = par_.ds->GetField(0);
        auto k = par_.ds->GetField(1);
        auto a = par_.ds->GetField(2);
        auto nx = par_.ds->GetField(3);
        auto ny = par_.ds->GetField(4);
        auto nz = par_.ds->GetField(5);

        using MIdx = typename M::MIdx;
        MIdx gs = bc.GetDimensions(); // global mesh size
        Scal ext = var.Double["extent"];
        Rect<Vect> d(Vect(0), Vect(gs) * (ext / gs.norminf())); // domain
        MIdx o(0); // cell origin

        auto gm = InitUniformMesh<M>(d, o, gs, 0);

        std::vector<std::vector<Vect>> vv;
        Vect h = GetCellSize(gm);
        for (auto c : gm.Cells()) {
          if (0.5 - std::abs(u[c] - 0.5) > 0.01) {
            auto ee = solver::GetLineEnds(Vect(nx[c], ny[c], nz[c]), a[c], h);
            Vect ea = ee[0];
            Vect eb = ee[1];
            auto xc = gm.GetCenter(c);
            Vect dz(0.,0.,h[2]*0.5);
            vv.push_back({xc+ea-dz, xc+eb-dz, xc+eb+dz, xc+ea+dz});
          }
        }
        auto fn = "s" + std::to_string(dmf_.GetN()) + ".vtk";
        WriteVtkPoly(vv, fn);
      }
    }
  }
}

template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("advection");
  Init(sem);
  Dump(sem);
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
      ff_flux_[f] = par_.fv(x, t).dot(m.GetSurface(f));
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
  if (sem("stat-loc")) {
    sumu_ = 0.;
    auto& u = as_->GetField();
    for (auto c : m.Cells()) {
      sumu_ += u[c] * m.GetVolume(c);
    }
    m.Reduce(&sumu_, "sum");
  }
  if (sem("stat")) {
    if (IsLead()) {
      ++var.Int["iter"];
      var.Double["t"] = as_->GetTime();
    }
    if (IsRoot()) {
      Scal t = var.Double["t"];
      Scal dt = var.Double["dt"];
      if (dms_.Try(t, dt)) {
        std::cout 
            << "t=" << t 
            << " dt=" << dt 
            << std::setprecision(16) << " sumu=" << sumu_
            << std::endl;
      }
    }
  }
  Dump(sem);
  sem.LoopEnd();
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  using K = Advection<M>;
  using Par = typename K::Par;
  Par par;
  par.fu0 = CreateInitU<M>(var);
  par.fv = CreateInitVel<Vect>(var);

  DistrSolver<M, K> ds(comm, var, par);
  par.ds = &ds;
  ds.Run();
}

int main(int argc, const char** argv) {
  //using Scal = double;
  //using Vect = GVect<Scal, 3>;
  //std::vector<std::vector<Vect>> vv{
  //  {Vect(0.,0.,0.), Vect(1.,0.,0.), Vect(1.,1.,0.), Vect(0.,1.,0.)}
  //};
  //WriteVtkPoly(vv, "s.vtk", "comment");

  return RunMpi(argc, argv, Main);
}
