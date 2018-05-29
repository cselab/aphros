#include <cassert>
#include <vector>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <iomanip>

#include "linear/hypre.h"

using MIdx = typename Hypre::MIdx;
using Scal = typename Hypre::Scal;
using Block = typename Hypre::Block;


bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-10;
}

// Print CMP if false
#define PFCMP(a, b) \
  if (!Cmp(a, b)) { \
    std::cerr \
      << std::scientific << std::setprecision(16) \
      << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
    assert(false);\
  }

int main (int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  size_t dim = 3;

  std::vector<MIdx> st = {{0, 0, 0}, {1, 0, 0}};

  int bs = 16;

  MIdx gs = {bs, bs, bs}; // global size

  Block b;
  b.l = {0, 0, 0};
  b.u = {bs-1, bs-1, bs-1};
  b.st = st;

  size_t n = 1;
  for (int i = 0; i < dim; ++i) {
    n *= b.u[i] - b.l[i] + 1;
  }

  // exact solution
  auto f = [](Scal x, Scal y, Scal z) {
    return std::sin(x);
  };

  std::vector<Scal> da(n * b.st.size());
  {
    size_t j = 0;
    for (int z = b.l[2]; z <= b.u[2]; ++z) {
      for (int y = b.l[1]; y <= b.u[1]; ++y) {
        for (int x = b.l[0]; x <= b.u[0]; ++x) {
          da[j] = 1.;
          ++j;
          da[j] = 0.5;
          ++j;
        }
      }
    }
    assert(j == n * b.st.size());
  }
  std::vector<Scal> dr(n);
  std::vector<Scal> dx(n);
  {
    size_t j = 0;
    for (int z = b.l[2]; z <= b.u[2]; ++z) {
      for (int y = b.l[1]; y <= b.u[1]; ++y) {
        for (int x = b.l[0]; x <= b.u[0]; ++x) {
          int xp = (x + 1 + gs[0]) % gs[0];
          dr[j] = da[2*j] * f(x, y, z) + da[2*j + 1] * f(xp, y, z);
          dx[j] = 0.;
          ++j;
        }
      }
    }
    assert(j == n);
  }

  b.a = &da;
  b.r = &dr;
  b.x = &dx;

  std::vector<Block> bb;
  bb.push_back(b);

  MPI_Comm comm = MPI_COMM_WORLD;
  std::vector<bool> per = {1, 0, 0};
  Scal tol = 1e-12;
  int print = 2;
  int maxiter = 80;

  Hypre h(comm, bb, gs, per, tol, print, "gmres", maxiter);
  h.Solve();

  // Check solution
  {
    size_t j = 0;
    for (int z = b.l[2]; z <= b.u[2]; ++z) {
      for (int y = b.l[1]; y <= b.u[1]; ++y) {
        for (int x = b.l[0]; x <= b.u[0]; ++x) {
          Scal e = f(x, y, z);
          PFCMP((*b.x)[j], e);
          ++j;
        }
      }
    }
  }

  MPI_Finalize();	
}
