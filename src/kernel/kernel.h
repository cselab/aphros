#pragma once

#include <string>
#include <array>

#include "parse/vars.h"

// TODO: remove h_gridpoint from MyBlockInfo
struct MyBlockInfo {
  using Idx = std::array<int, 3>;
  Idx index;
  void* ptrBlock;
  double h_gridpoint;
  double origin[3];
  Idx bs;
  int hl; // number of halo cells
  bool isroot; // root block (one among blocks on all PEs)
  bool islead; // lead block (one per each PE)
  Idx gs; // global size
  size_t maxcomm; // maximum number of communication requests
};

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

