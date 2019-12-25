#undef NDEBUG
#include <mpi.h>
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "distr/distrsolver.h"
#include "kernel/kernelmeshpar.h"
#include "parse/vars.h"

#include "dump/dumper.h"
#include "dump/output.h"
#include "func/init_u.h"
#include "func/init_vel.h"
#include "geom/mesh.h"
#include "geom/vect.h"
#include "parse/tvd.h"
#include "parse/vof.h"
#include "solver/advection.h"
#include "solver/solver.h"
#include "solver/tvd.h"
#include "solver/vof.h"
#include "solver/vofm.h"
#include "util/suspender.h"

#include "dump/dump.h"
#include "dump/vtk.h"

#include "cmp.h"

using namespace solver;

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

  std::function<Vect(Vect, Scal)> fv; // velocity

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
  using ASVM = solver::Vofm<M>; // advection VOF

  // using P::P; // inherit constructor
  Advection(Vars& var, const MyBlockInfo& bi, Par& par);
  void Run() override;
  void Init(Sem& sem);
  void Dump(Sem& sem);

 protected:
  using P::IsLead;
  using P::IsRoot;
  using P::m;
  using P::par_;
  using P::var;

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

  // boundary conditions for advection (empty)
  MapCondFaceAdvection<Scal> bc_;
};

template <class M>
Advection<M>::Advection(Vars& v, const MyBlockInfo& b, Par& p)
    : KernelMeshPar<M, Par>(v, b, p)
    , dmf_(v, "dump_field_")
    , dms_(v, "dump_stat_") {}

template <class M>
void Advection<M>::Init(Sem& sem) {
  if (sem.Nested("init-field")) {
    InitVf(fcu_, var, m);
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

    // cell conditions for advection (empty)
    MapCell<std::shared_ptr<solver::CondCell>> mc_cond;

    // source
    fc_src_.Reinit(m, 0.);

    std::string as = var.String["advection_solver"];
    if (as == "tvd") {
      auto p = std::make_shared<typename AST::Par>();
      Parse<M>(p.get(), var);
      as_.reset(
          new AST(m, fcu_, bc_, &ff_flux_, &fc_src_, 0., var.Double["dt"], p));
    } else if (as == "vof") {
      auto p = std::make_shared<typename ASV::Par>();
      Parse<M, ASV>(p.get(), var);
      p->dmp = std::unique_ptr<Dumper>(new Dumper(var, "dump_part_"));
      auto fccl = FieldCell<Scal>(m, 0);
      as_.reset(new ASV(
          m, fcu_, fccl, bc_, &ff_flux_, &fc_src_, 0., var.Double["dt"], p));
    } else if (as == "vofm") {
      auto p = std::make_shared<typename ASVM::Par>();
      Parse<M, ASVM>(p.get(), var);
      p->dmp = std::unique_ptr<Dumper>(new Dumper(var, "dump_part_"));
      auto fccl = FieldCell<Scal>(m, 0);
      as_.reset(new ASVM(
          m, fcu_, fccl, bc_, &ff_flux_, &fc_src_, 0., var.Double["dt"], p));
    } else {
      throw std::runtime_error("Unknown advection_solver=" + as);
    }
  }
}

template <class M, class Vect = typename M::Vect>
Vect GetCellSize(const M& m) {
  Vect h; // result
  IdxCell c0(0);
  h = m.GetNode(m.GetNode(c0, 7)) - m.GetNode(m.GetNode(c0, 0));
  assert(std::abs(h.prod() - m.GetVolume(c0)) < 1e-10);
  return h;
}

template <class M>
void Advection<M>::Dump(Sem& sem) {
  // TODO: Suspender: allow change stages between time steps
  if (sem("dump")) {
    if (dmf_.Try(var.Double["t"], var.Double["dt"])) {
      m.Dump(&as_->GetField(), "u");
      if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
        m.Dump(&as->GetAlpha(), "a");
        auto& n = as->GetNormal();
        m.Dump(&n, 0, "nx");
        m.Dump(&n, 1, "ny");
        m.Dump(&n, 2, "nz");
        m.Dump(&as->GetCurv(), "k");
      }
      if (auto as = dynamic_cast<solver::Vofm<M>*>(as_.get())) {
        m.Dump(as->GetAlpha()[0], "a");
        auto* n = as->GetNormal()[0];
        m.Dump(n, 0, "nx");
        m.Dump(n, 1, "ny");
        m.Dump(n, 2, "nz");
        m.Dump(&as->GetCurvSum(), "k");
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
}

template <class M>
void Advection<M>::Run() {
  auto sem = m.GetSem("advection");
  Init(sem);
  Dump(sem);
  sem.LoopBegin();
  if (sem("empty")) {
    // nop // TODO: bugfix loop, empty stage needed
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
        auto f = m.GetFace(c, q);
        maxvel_ = std::max(maxvel_, std::abs(ff_flux_[f]) / m.GetVolume(c));
      }
    }
    m.Reduce(&maxvel_, "max");
  }
  if (sem("cfl")) {
    Scal dt = var.Double["cfl"] / maxvel_;
    dt = std::min(dt, var.Double["dtmax"]);
    as_->SetTimeStep(dt);
    if (m.IsLead()) {
      this->var_mutable.Double.Set("dt", as_->GetTimeStep());
    }
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
    if (m.IsLead()) {
      ++(this->var_mutable.Int["iter"]);
      this->var_mutable.Double.Set("t", as_->GetTime());
    }
    if (IsRoot()) {
      Scal t = var.Double["t"];
      Scal dt = var.Double["dt"];
      if (dms_.Try(t, dt)) {
        std::cout << "t=" << t << " dt=" << dt << std::setprecision(16)
                  << " sumu=" << sumu_ << std::endl;
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
  par.fv = CreateInitVel<Vect>(var);

  DistrSolver<M, K> ds(comm, var, par);
  par.ds = &ds;
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
