#undef NDEBUG
#include <mpi.h>
#include <fstream>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using MIdx = typename M::MIdx;

struct State {};

std::ofstream out;

void Run(M& m, State&, Vars&) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fc;
  }* ctx(sem);
  auto& fc = ctx->fc;
  auto& bc = m.GetInBlockCells();
  MIdx w = bc.GetBegin() / bc.GetSize();
  if (sem()) {
    if (m.IsLead()) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      out.open(std::to_string(rank) + ".out");
    }
  }
  if (sem("A")) {
    fc.Reinit(m);
    out << "A" << w << std::endl;
    m.Comm(&fc);
  }
  if (sem("B")) {
    out << "B" << w << std::endl;
  }
  // FIXME: without empty stage:
  // aborted _ZNK15SynchronizerMPIIN13cubismnc_impl13GFieldViewRawINS0_
  if (sem()) {}
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 3
set int by 1
set int bz 3

set int bsx 8
set int bsy 8
set int bsz 8

set int px 2
set int py 1
set int pz 2

set int openmp 0
set int histogram 0
set int mpi_compress_msg 0
)EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

