#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>

#include "geom/mesh.h"
#include "solver/solver.h"
#include "func/init_u.h"
#include "parse/vars.h"
#include "solver/partstrmesh.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"

using namespace solver;

const int dim = 3;
using MIdx = GMIdx<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;

Mesh GetMesh(MIdx s /*size in cells*/) {
  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx b(0, 0, 0); // lower index
  int hl = 2;         // halos 
  return InitUniformMesh<Mesh>(dom, b, s, hl, true, true, s, 0);
}

Mesh* _mesh;

struct GPar {};

template <class M_>
class Simple : public KernelMeshPar<M_, GPar> {
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
  using P::bi_;
  using P::m;

 private:
  FieldCell<Scal> fcu_;
  FieldCell<Scal> fca_;
  FieldCell<Scal> fck_;
  FieldCell<Vect> fcn_;
  FieldCell<bool> fci_;
  std::shared_ptr<solver::PartStrMesh<M>> psm_;
};


#include "chpartstr.h"

template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem();
  if (sem()) {
    _mesh = &m;

    fcu_.Reinit(m);
    fca_.Reinit(m);
    fcn_.Reinit(m);

    Vars par;
    par.SetStr("string", "init_vf", "list");
    par.SetStr("string", "list_path", "b.dat");
    par.SetStr("int", "list_ls", "2");
    par.SetStr("int", "dim", "3");

    auto fu = CreateInitU<Mesh>(par);
    fu(fcu_, m);

    std::cout << GetNcInter(fcu_);

    FieldCell<Scal> nx(m), ny(m), nz(m);
    vector nn = {nx, ny, nz};
    CalcNormal(fcu_, nn);

    DumpFacets(fcu_, "o.vtk");
    DumpLines(fcu_, nn, kPartstr, "line.vtk");

    kPartstr.csv = true;

    for (auto c : m.Cells()) {
      if (fcu_[c] > 0 && fcu_[c] < 1) {
        std::cout << partstr(c, fcu_, nn) << std::endl;
      }
    }

    auto& kp = kPartstr;
    auto ps = std::make_shared<typename PartStr<Scal>::Par>();
    ps->npmax = kp.Np;
    ps->relax = kp.eta;
    ps->leq = kp.Hp;

    auto pm = std::make_shared<typename solver::PartStrMesh<Mesh>::Par>();
    pm->tol = kp.eps;
    pm->ns = kp.Ns;
    pm->itermax = kp.itermax;
    pm->ps = ps;

    psm_ = std::make_shared<solver::PartStrMesh<Mesh>>(m, pm);

    fci_ = DetectInterface(fcu_);
    solver::UNormal<Mesh>::CalcNormal(m, fcu_, fci_, dim, fcn_);
    for (auto c : m.SuCells()) {
      if (fci_[c]) {
        fca_[c] = Reconst<Scal>::GetLineA(fcn_[c], fcu_[c], m.GetCellSize());
      }
    }
    fck_.Reinit(m, GetNan<Scal>());
  }

  if (sem.Nested()) {
    psm_->Part(fcu_, fca_, fcn_, fci_, fck_,
             MapFace<std::shared_ptr<solver::CondFace>>());
  }

  if (sem.Nested()) {
    psm_->DumpParticles(fca_, fcn_, -1, 0);
  }
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using K = Simple<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}


int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
