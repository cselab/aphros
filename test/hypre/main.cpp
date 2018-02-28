#include <cassert>
#include <vector>
#include <iostream>
#include <mpi.h>

#include "hypre.h"


int main (int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  using MIdx = typename Hypre::MIdx;
  using Scal = typename Hypre::Scal;
  using Block = typename Hypre::Block;

  size_t dim = 3;

  std::vector<MIdx> st = {{0, 0, 0}, {-1, 0, 0}};

  int bs = 16;

  Block b;
  b.l = {0, 0, 0};
  b.u = {bs-1, bs-1, bs-1};
  b.st = st;

  size_t n = 1;
  for (int i = 0; i < dim; ++i) {
    n *= b.u[i] - b.l[i] + 1;
  }

  std::vector<Scal> a(n * b.st.size());
  std::vector<Scal> r(n);
  std::vector<Scal> x(n);

  b.a = &a;
  b.r = &r;
  b.x = &x;

  std::vector<Block> bb;
  bb.push_back(b);

  MPI_Comm comm = MPI_COMM_WORLD;
  MIdx gs = {bs, bs, bs};
  std::vector<bool> per = {0, 0, 0};
  Scal tol = 1e-10;
  int print = 2;

  Hypre h(comm, bb, gs, per, tol, print);
  h.Solve();

  MPI_Finalize();	
}
