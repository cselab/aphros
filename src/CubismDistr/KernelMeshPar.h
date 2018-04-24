#pragma once

#include <memory>

#include "hydro/vect.hpp"
#include "hydro/linear.hpp"
#include "Kernel.h"
#include "Vars.h"
#include "KernelMesh.h"

// Factory for ParKernel
// M_: mesh
// K_: kernel derived from KernelMesh
template <class M_, class K_>
class KernelMeshParFactory : public KernelMeshFactory<M_> {
 public:
  using M = M_;
  using K = K_;
  GKernelMeshFactory(std::function<Scal(Vect)> fu0,
                   std::function<Vect(Vect,Scal)> fvel) 
      : fu0_(fu0), fvel_(fvel) {}
  K* Make(Vars& par, const MyBlockInfo& bi) override {
    return new K(par, bi, fu0_, fvel_);
  }
};

