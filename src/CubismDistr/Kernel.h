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
  bool isroot; // root block (single among all PEs)
  bool islead; // lead block (one per each PE)
};

// Suspendable kernel
class Kernel {
 public:
  virtual void Run() = 0;
  virtual ~Kernel() {}
};

class KernelFactory {
 public:
  virtual ~KernelFactory() {}
  virtual Kernel* Make(Vars&, const MyBlockInfo&) = 0;
};

