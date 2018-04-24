#pragma once

#include <memory>
#include <mpi.h>
#include <stdexcept>

#include "CubismDistr/Vars.h"
#include "CubismDistr/KernelMeshPar.h"
#include "CubismDistr/Distr.h"
#include "CubismDistr/ICubism.h"
#include "CubismDistr/ILocal.h"
#include "CubismDistr/Interp.h"

#include "hydro/vect.hpp"
#include "hydro/block.h"

using namespace geom;

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
    Distr* b; // base
    if (var.Int["loc"]) {
      b = CreateLocal(comm, kf, var);
    } else {
      b = CreateCubism(comm, kf, var);
    }

    d_.reset(dynamic_cast<D*>(b));
    if (!d_) {
      throw std::runtime_error("DistrSolver: Can't cast to D");
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
  FieldCell<Scal> GetField() const {
    return d_->GetGlobalField(0);
  }
  geom::FieldCell<Scal> GetField(size_t i) const {
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
           std::function<void(MPI_Comm, Vars&)> r) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool isroot = (!rank);

  Vars var;   // parameter storage
  Interp ip(var); // interpretor (parser)

  std::string fn = "a.conf";
  if (argc == 1) {
    // nop
  } else if (argc == 2) {
    fn = argv[1];
  } else {
    if (isroot) {
      std::cerr << "usage: " << argv[0] << " [a.conf]" << std::endl;
    }
    return 1;
  }

  if (isroot) {
    std::cerr << "Loading config from '" << fn << "'" << std::endl;
  }

  std::ifstream f(fn);  // config file
  // Read file and run all commands
  ip.RunAll(f);   

  // Print vars on root
  if (isroot) {            
    std::cerr << "\n=== config begin ===" << std::endl;
    ip.PrintAll(std::cerr);
    std::cerr << "=== config end ===\n" << std::endl;
  }

  bool loc = var.Int["loc"];

  MPI_Comm comm;
  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      r(comm, var);
    }
  } else {
    comm = MPI_COMM_WORLD;
    r(comm, var);
  }

  MPI_Finalize();	
  return 0;
}
