// Created by Petr Karnakov on 27.11.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "linear/hypre.h"
#include "linear/hypresub.h"
#include "util/histogram.h"

using MIdx = typename Hypre::MIdx;
using Scal = typename Hypre::Scal;
using Block = typename Hypre::Block;

#define EV(x) (#x) << "=" << (x) << " "

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-10;
}

// Print CMP if false
#define PFCMP(a, b)                                                         \
  if (!Cmp(a, b)) {                                                         \
    std::cerr << std::scientific << std::setprecision(16) << #a << "=" << a \
              << ", " << #b << "=" << b << std::endl;                       \
    assert(false);                                                          \
  }

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  constexpr size_t dim = 3;

  std::vector<MIdx> stencil = {{0, 0, 0}, {1, 0, 0}};

  int bs = 16;

  MIdx gs = {bs, bs, bs}; // global size

  Block b;
  b.l = {0, 0, 0};
  b.u = {bs - 1, bs - 1, bs - 1};
  b.stencil = stencil;

  size_t n = 1;
  for (size_t i = 0; i < dim; ++i) {
    n *= b.u[i] - b.l[i] + 1;
  }

  // exact solution
  auto f = [](Scal x, Scal y, Scal z) {
    return std::sin(x) * std::cos(y) * std::exp(z / 100.);
  };

  std::vector<Scal> da(n * b.stencil.size());
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
    assert(j == n * b.stencil.size());
  }
  std::vector<Scal> dr(n);
  std::vector<Scal> dx(n);
  {
    size_t j = 0;
    for (int z = b.l[2]; z <= b.u[2]; ++z) {
      for (int y = b.l[1]; y <= b.u[1]; ++y) {
        for (int x = b.l[0]; x <= b.u[0]; ++x) {
          int xp = (x + 1 + gs[0]) % gs[0];
          dr[j] = da[2 * j] * f(x, y, z) + da[2 * j + 1] * f(xp, y, z);
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
  bb.push_back(b);
  bb.push_back(b);

  MPI_Comm comm = MPI_COMM_WORLD;
  MIdx per = {1, 0, 0};
  Scal tol = 1e-12;
  int print = 0;
  int maxiter = 80;

  (void)tol;
  (void)print;
  (void)maxiter;

  int rank;
  MPI_Comm_rank(comm, &rank);

  MPI_Comm commsub;
  MPI_Comm_split(comm, rank, rank, &commsub);

  int ranksub;
  MPI_Comm_rank(commsub, &ranksub);
  std::cout << EV(rank) << EV(ranksub) << std::endl;

  HypreSub::InitServer(comm, commsub);

  {
    Histogram hist(comm, "hypresub", true);
    if (ranksub == 0) {
      {
        HypreSub d(comm, bb, gs, per);
        d.Solve(tol, print, "gmres", maxiter);
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
      }
      HypreSub::StopServer();
    } else {
      HypreSub::RunServer(hist);
    }
  }

  MPI_Finalize();
}
