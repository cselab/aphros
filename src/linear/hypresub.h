// Created by Petr Karnakov on 27.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "hypre.h"
#include "util/histogram.h"

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
  HypreSub(MPI_Comm comm, const std::vector<Block>& bb, MIdx gs, MIdx per);
  HypreSub() = delete;
  HypreSub(const HypreSub&) = delete;
  ~HypreSub();
  // Initializes server
  // comm: full communicator
  // commsub: communicator over ranks on the same node
  // Needs to be called by all ranks
  static void InitServer(MPI_Comm comm, MPI_Comm commsub);
  // Runs server loop.
  // Needs to be called by non-zero ranks in commsub
  static void RunServer(Sampler& samp);
  static void Send(std::string cmd, int rank);
  static void Send(const Block& b, int rank);
  static void Send(const std::vector<Block>& bb, int rank);
  // Stops server on non-zero ranks in commsub
  // Needs to be called by zero ranks in commsub
  static void StopServer();

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
