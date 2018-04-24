#pragma once

#include <vector>
#include <cassert>
#include <mpi.h>
#include <memory>
#include <string>

// TODO: Convention of *_ for private variables ignored

class Hypre {
 public:
  using MIdx = std::vector<int>; 
  using Scal = double;
  struct Block { // linear system ax=r
    MIdx l; // lower corner
    MIdx u; // upper corner
    std::vector<MIdx> st; // stencil
    std::vector<Scal>* a; // matrix coeffs of size n * st.size()
    std::vector<Scal>* r; // rhs of size n
    std::vector<Scal>* x; // solution and initial guess of size n
  };

  Hypre(MPI_Comm comm, 
        std::vector<Block> bb, 
        MIdx gs /*global size*/,
        std::vector<bool> per /*periodic per direction*/,
        Scal tol /*tolerance*/,
        int print /*print level*/,
        std::string solver,
        int maxiter);
  Hypre() = delete;
  Hypre(const Hypre&) = delete;

  void Solve();
  ~Hypre();
 
 private:
  size_t dim;
  std::vector<Block> bb;
  struct HypreData;
  std::unique_ptr<HypreData> hd;
  std::string solver_;
  int maxiter_;
};
