#pragma once

#include <memory>
#include <mpi.h>
#include <stdexcept>
#include <functional>

#include "parse/vars.h"
#include "kernel/kernelmeshpar.h"
#include "distr.h"
#include "cubism.h"
#include "cubismnc.h"
#include "local.h"
#include "parse/parser.h"
#include "geom/vect.h"
#include "geom/block.h"

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
  using D = DistrMesh<KernelMeshFactory<M>>;
  using Par = typename K::Par;

  DistrSolver(MPI_Comm comm, Vars& var, Par& par)
      : var(var), par_(par)
  {
    // Create kernel factory
    KF kf(par_);

    var.Double.Set("t", 0.);

    // Initialize blocks
    std::unique_ptr<Distr> b;

    const std::string be = var.String["backend"];
    std::cerr << "backend=" << be << std::endl;
    if (be == "local") {
      std::cerr << "CreateLocal" << std::endl;
      b = CreateLocal(comm, kf, var);
    } else if (be == "cubism") {
      std::cerr << "CreateCubism" << std::endl;
      b = CreateCubism(comm, kf, var);
    } else if (be == "cubismnc") {
      std::cerr << "CreateCubismnc" << std::endl;
      b = CreateCubismnc(comm, kf, var);
    } else {
      throw std::runtime_error("DistrSolver: unknown backend='" + be + "'");
    }

    d_ = std::unique_ptr<D>(static_cast<D*>(b.release()));
    if (!d_) {
      throw std::runtime_error("DistrSolver: no compatible backend found");
    }
  }
  void Run() {
    d_->Run();
  }
  void Report() {
    d_->Report();
  }
  GBlock<IdxCell, dim> GetBlock() const {
    return d_->GetGlobalBlock();
  }
  GIndex<IdxCell, dim> GetIndex() const {
    return d_->GetGlobalIndex();
  }
  FieldCell<Scal> GetField() const {
    return d_->GetGlobalField(0);
  }
  FieldCell<Scal> GetField(size_t i) const {
    return d_->GetGlobalField(i);
  }
  double GetTime() const { 
    return var.Double["t"]; 
  }
  double GetTimeStep() const { 
    return var.Double["dt"]; 
  }

 private:
  Vars& var;
  std::unique_ptr<D> d_;
  Par& par_;
};

int RunMpi(int argc, const char ** argv,
           std::function<void(MPI_Comm, Vars&)> r);
