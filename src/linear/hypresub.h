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
  HypreSub(MPI_Comm comm, const std::vector<Block>& bb, MIdx gs, MIdx per);
  HypreSub() = delete;
  HypreSub(const HypreSub&) = delete;
  ~HypreSub();
  // Initializes server
  // comm: full communicator
  // commsub: sub-communicator (one rank in each node)
  // Needs to be called by ranks in commsub.
  static void InitServer(MPI_Comm comm, MPI_Comm commsub);
  // Runs server loop.
  // Needs to be called by all ranks in comm which are not in commsub.
  static void RunServer();
  static void Send(std::string cmd, int rank);
  static void Send(const Block& b, int rank);
  static void Send(const std::vector<Block>& bb, int rank);
  // Stops server on all ranks which are not in commsub
  // Needs to be called by ranks in commsub.
  static void StopServer();

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
