#pragma once

#include <string>
#include <array>

#include "hydro/suspender.h"
#include "Vars.h"

using Scal = double;

struct MyElem {
  static const size_t s = 8;
  Scal a[s];
  void init(Scal val) { 
    for (size_t i = 0; i < s; ++i) {
      a[i] = val;
    }
  }
  void clear() {
    init(0);
  }

  MyElem& operator=(const MyElem&) = default;
};


struct MyBlock {
  static const int bs = _BLOCKSIZE_;
  static const int sx = bs;
  static const int sy = bs;
  static const int sz = bs;
  static const int n = sx * sy * sz;

  // floats per element
  static const int fe = sizeof(MyElem) / sizeof(Scal);
};


// TODO: remove h_gridpoint from MyBlockInfo
struct MyBlockInfo {
  using Idx = std::array<int, 3>;
  Idx index;
  void* ptrBlock;
  double h_gridpoint;
  double origin[3];
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

class KernelFactory {
 public:
  virtual Kernel* Make(Vars&, const MyBlockInfo&) = 0;
  virtual ~KernelFactory() {}
};





