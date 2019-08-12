#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>

#include "geom/mesh.h"
#include "solver/solver.h"
#include "parse/vars.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
using namespace solver;


const int dim = 3;
using MIdx = GMIdx<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;

extern "C" {
#include "imp.h"
TFunc kFunc;
}

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
};



template <class M>
void Simple<M>::Run() {
  auto sem = m.GetSem();
  if (sem()) {
    fcu_.Reinit(m);
    m.Dump(&fcu_, "u");
  }
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using K = Simple<M>;
  using Par = typename K::Par;
  Par par;

  double fcu[] = {1, 2, 3};
  kFunc(var.Int["a"], var.Int["b"], fcu);
  return;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}


int CMain(TFunc f) {
  int argc = 1;
  const char* argv[] = {"main"};
  kFunc = f;
  return RunMpi(argc, argv, Main);
}

