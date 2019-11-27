#pragma once

#include <vector>
#include <array>
#include <mpi.h>
#include <memory>
#include <string>

#include "hypre.h"

class HypreSub {
 public:
  static constexpr size_t dim = Hypre::dim;
  using MIdx = Hypre::MIdx;
  using Scal = Hypre::Scal;
  using Block = Hypre::Block;

  // bb: blocks
  // gs: global size
  // per: periodic in each direction
  // tol: tolerance
  // print: print level
  // solver: solver name
  // maxiter: maximum number of iterations
  HypreSub(MPI_Comm comm, const std::vector<Block>& bb,
           MIdx gs, std::vector<bool> per);
  HypreSub() = delete;
  HypreSub(const HypreSub&) = delete;
  ~HypreSub();

  // Assembles matrix and vectors from bb
  void Update();
  // Solves system and puts result to x
  void Solve(Scal tol, int print, std::string solver, int maxiter);
  // Returns relative residual norm from last Solve()
  Scal GetResidual();
  // Returns the number of iterations from last Solve()
  int GetIter();

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
