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
#include "parse/vof.h"
#include "solver/tvd.h"
#include "parse/tvd.h"
#include "func/init_vel.h"
#include "func/init_u.h"
#include "dump/dumper.h"
#include "dump/output.h"

#include "dump/dump.h"
#include "dump/vtk.h"

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
  using AST = solver::Tvd<M>; // advection TVD
  using ASV = solver::Vof<M>; // advection VOF

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
    int vdim = var.Int["dim"];
    for (auto f : m.AllFaces()) {
      Vect x = m.GetCenter(f);
      if (vdim == 2) {
        x[2] = 0.;
      }
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
      auto p = std::make_shared<typename AST::Par>();
      Parse<M>(p.get(), var);
      as_.reset(new AST(m, fcu_, bc, &ff_flux_, 
                       &fc_src_, 0., var.Double["dt"], p));
    } else if (as == "vof") {
      auto p = std::make_shared<typename ASV::Par>();
      Parse<M>(p.get(), var);
      p->dmp = std::unique_ptr<Dumper>(new Dumper(var, "dump_part_"));
      as_.reset(new ASV(m, fcu_, bc, &ff_flux_, 
                       &fc_src_, 0., var.Double["dt"], p));
    } else {
      throw std::runtime_error("Unknown advection_solver=" + as);
    }
  }
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
  // TODO: Suspender: allow change stages between time steps
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

        auto& kh = const_cast<FieldCell<Scal>&>(as->GetCurvH());
        m.Dump(&kh, "kh");
        auto& kp = const_cast<FieldCell<Scal>&>(as->GetCurvP());
        m.Dump(&kp, "kp");
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
      if (dmf_.Try(var.Double["t"], var.Double["dt"])) {
        if (IsLead()) {
          auto u = par_.ds->GetField(0);
          auto k = par_.ds->GetField(1);
          auto a = par_.ds->GetField(2);
          auto nx = par_.ds->GetField(3);
          auto ny = par_.ds->GetField(4);
          auto nz = par_.ds->GetField(5);
        }
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
  if (sem.Nested("post")) {
    as_->PostStep();
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
  return RunMpi(argc, argv, Main);
}
