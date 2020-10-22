// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
#include <functional>
#include <memory>
#include <stdexcept>

#include "cubismnc.h"
#include "distr.h"
#include "geom/block.h"
#include "geom/vect.h"
#include "kernel/kernelmeshpar.h"
#include "local.h"
#include "parse/parser.h"
#include "parse/vars.h"

// Client for DistrMesh
// M_: mesh
// K_: kernel derived from KernelMeshPar
template <class M_, class K_>
class DistrSolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  using K = K_;
  using KF = KernelMeshParFactory<M, K>;
  using Par = typename K::Par;

  DistrSolver(MPI_Comm comm, Vars& var0, Par& par)
      : var(var0), var_mutable(var0), kernelfactory_(par) {
    const std::string be = var.String["backend"];
    if (be == "local") {
      d_ = CreateLocal<M>(comm, kernelfactory_, var_mutable);
    } else if (be == "cubismnc") {
      d_ = CreateCubismnc<M>(comm, kernelfactory_, var_mutable);
    } else {
      throw std::runtime_error("DistrSolver: unknown backend='" + be + "'");
    }
  }
  void Run() {
    d_->Run();
  }
  void Report() {
    d_->Report();
  }

 private:
  const Vars& var;
  Vars& var_mutable;
  std::unique_ptr<DistrMesh<M>> d_;
  KF kernelfactory_;
};

class MpiWrapper {
 public:
  MpiWrapper(int* argc, const char*** argv);
  ~MpiWrapper();
  MPI_Comm GetComm() const;
  static int GetCommSize(MPI_Comm);
  static int GetCommRank(MPI_Comm);
  static bool IsRoot(MPI_Comm);
  int GetCommSize() const;
  int GetCommRank() const;
  bool IsRoot() const;

 private:
  MPI_Comm comm_;
};

#define MPICALL(x)                                                    \
  do {                                                                \
    int errorcode;                                                    \
    errorcode = x;                                                    \
    if (errorcode != MPI_SUCCESS) {                                   \
      char string[MPI_MAX_ERROR_STRING];                              \
      int resultlen;                                                  \
      MPI_Error_string(errorcode, string, &resultlen);                \
      throw std::runtime_error(FILELINE + ": mpi failed: " + string); \
    }                                                                 \
  } while (0)

int RunMpi(int argc, const char** argv, std::function<void(MPI_Comm, Vars&)> r);
