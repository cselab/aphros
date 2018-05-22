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
  using Sem = typename Mesh::Sem;

  //using P::P; // inherit constructor
  Advection(Vars& var, const MyBlockInfo& bi, Par& par);
  void Run() override;
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
  std::unique_ptr<solver::AdvectionSolver<M>> as_; // advection solver
  std::function<Vect(Vect,Scal)> fvel_;
  Scal maxvel_; // maximum velocity relative to cell length [1/time]
                // cfl = dt * maxvel
  Scal sumu_; // sum of fluid volume
  FieldCell<Scal> fcnx_, fcny_, fcnz_; // normal to interface (tmp)
                                             // used for Vof dump
  Dumper dumper_;
};

template <class M>
Advection<M>::Advection(Vars& var, const MyBlockInfo& bi, Par& par)
    : KernelMeshPar<M, Par>(var, bi, par)
    , dumper_(var)
{
  // initial field for advection
  FieldCell<Scal> u0(m, 0);
  for (auto c : m.AllCells()) {
    Vect x = m.GetCenter(c);
    u0[c] = par_.fu0(x);
  }

  // boundary conditions for advection (empty)
  MapFace<std::shared_ptr<solver::CondFace>> bc;

  // cell conditions for advection (empty)
  MapCell<std::shared_ptr<solver::CondCell>> mc_cond;

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
    using AS = solver::AdvectionSolverExplicit<M>;
    auto p = std::make_shared<typename AS::Par>();
    p->sharp = var.Double["sharp"];
    p->sharpo = var.Double["sharpo"];
    p->split = var.Int["split"];
    as_.reset(new AS(m, u0, bc, &ff_flux_, 
                     &fc_src_, 0., var.Double["dt"], p));
  } else if (as == "vof") {
    using AS = solver::Vof<M>;
    auto p = std::make_shared<typename AS::Par>();
    as_.reset(new AS(m, u0, bc, &ff_flux_, 
                     &fc_src_, 0., var.Double["dt"], p));
  } else {
    throw std::runtime_error("Unknown advection_solver=" + as);
  }
}

template <class M>
void Advection<M>::Dump(Sem& sem) {
  if (sem("dump")) {
    if (dumper_.Try(var.Double["t"], var.Double["dt"])) {
      auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
      m.Dump(&u, "u");
      auto& k = const_cast<FieldCell<Scal>&>(as_->GetCurv());
      m.Dump(&k, "k");
      if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
        auto& a = const_cast<FieldCell<Scal>&>(as->GetAlpha());
        m.Dump(&a, "a");
        auto &n = as->GetNormal();
        fcnx_ = GetComponent(n, 0);
        fcny_ = GetComponent(n, 1);
        fcnz_ = GetComponent(n, 2);
        m.Dump(&fcnx_, "nx");
        m.Dump(&fcny_, "ny");
        m.Dump(&fcnz_, "nz");
      }

      if (IsRoot()) {
        dumper_.Report();
      }
    }
  }
  if (sem("dumpwrite")) {
    // Empty stage for DumpWrite
    // TODO: revise
  }
  if (sem("dump")) {
    if (dumper_.Try(var.Double["t"], var.Double["dt"])) {
      auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
      m.Dump(&u, "u");
      if (IsRoot()) {
        dumper_.Report();
      }
    }
  }
  if (sem("dumpwrite")) {
    // Empty stage for DumpWrite
    // TODO: revise
  }
}

template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("advection");
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
      sumu_ += u[c] * m.GetVolume(c);
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
  par.fu0 = CreateInitU<Vect>(var);
  par.fvel = CreateInitVel<Vect>(var);

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
