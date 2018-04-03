#pragma once

#include <string>
#include <array>

#include "hydro/suspender.h"
#include "Vars.h"


// TODO: remove h_gridpoint from MyBlockInfo
struct MyBlockInfo {
  using Idx = std::array<int, 3>;
  Idx index;
  void* ptrBlock;
  double h_gridpoint;
  double origin[3];
  Idx bs;
  int hl; // number of halo cells
};

// Suspendable kernel
class Kernel {
 public:
  virtual void Run() = 0;
  //virtual void ReadBuffer(LabMPI&) = 0;
  //virtual void WriteBuffer(Block_t&) = 0;
  virtual ~Kernel() {}
};

// TODO: revise isroot and islead
class KernelFactory {
 public:
  virtual Kernel* Make(Vars&, const MyBlockInfo&, 
                       bool isroot, bool islead) = 0;
  virtual ~KernelFactory() {}
};





