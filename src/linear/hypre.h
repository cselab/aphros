// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <HYPRE_struct_ls.h>

class Hypre {
 public:
  static constexpr size_t dim = 3;
  using MIdx = std::array<HYPRE_Int, dim>;
  using Scal = double;
  struct Block { // linear system ax=r
    MIdx l; // lower corner
    MIdx u; // upper corner
    std::vector<MIdx> stencil; // stencil
    std::vector<Scal>* a; // matrix coeffs of size n * stencil.size()
    std::vector<Scal>* r; // rhs of size n
    std::vector<Scal>* x; // solution and initial guess of size n
  };

  // bb: blocks
  // gs: global size
  // per: periodic in each direction
  // tol: tolerance
  // print: print level
  // solver: solver name
  // maxiter: maximum number of iterations
  Hypre(MPI_Comm comm, const std::vector<Block>& bb, MIdx gs, MIdx per);
  Hypre() = delete;
  Hypre(const Hypre&) = delete;
  ~Hypre();

  // Assembles matrix and vectors from bb
  void Update();
  // Solves system and puts result to x
  void Solve(Scal tol, int print, std::string solver, int maxiter);
  // Returns relative residual norm from last Solve()
  Scal GetResidual() const;
  // Returns the number of iterations from last Solve()
  int GetIter() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
