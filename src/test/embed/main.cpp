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
#include "solver/solver.h"
#include "solver/embed.h"
#include "dump/output.h"
#include "dump/output_paraview.h"

struct GPar {};

template <class M_>
class KernelEmbed : public KernelMeshPar<M_, GPar> {
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

 private:
  using EM = solver::Embed<M>;
  std::unique_ptr<EM> em_;
  FieldCell<Scal> fct_;  // cell type
};


template <class M>
void KernelEmbed<M>::Run() {
  auto sem = m.GetSem("Run");
  if (sem("init")) {
   FieldNode<Scal> fnf(m, 0);
   for (auto n : m.Nodes()) {
     auto x = m.GetNode(n);
     fnf[n] = x.norm() * (1 + std::sin(x[0] * 5)*0.1 + std::sin(x[1] * 7)*0.15) - 1.;
     fnf[n] *= -1;
   }
   em_ = std::unique_ptr<EM>(new EM(m, fnf));
   fct_.Reinit(m);
   for (auto c : m.Cells()) {
     fct_[c] = size_t(em_->GetCellType()[c]);
   }
   m.Dump(&fct_, "type");
  }
  if (sem.Nested("dump")) {
    em_->DumpPoly();
  }
}

using Scal = double;
using M = MeshStructured<Scal, 3>;
using K = KernelEmbed<M>;

void Main(MPI_Comm comm, Vars& var0) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  Vars var = var0;

  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}


int main (int argc, const char ** argv) {
  return RunMpi(argc, argv, Main);
}
