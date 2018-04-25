#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
#include "kernel/hydro.h"

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using K = Hydro<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}


int main (int argc, const char ** argv) {
  return RunMpi(argc, argv, Main);
}
