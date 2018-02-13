#pragma once

#include <array>
#include <mpi.h>
#include <memory>

#include "Kernel.h"

using Idx = std::array<int, 3>;

class Distr {
 public:
  virtual bool IsDone() const = 0;
  virtual void Step() = 0;
};

std::unique_ptr<Distr> CreateCubism(
    MPI_Comm comm, KernelFactory& kf, 
    int bs, Idx b, Idx p, int es, int h);
