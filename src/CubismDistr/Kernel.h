#pragma once

#include <string>
#include <array>

#include "hydro/suspender.h"
#include "Vars.h"

struct MyElem {
  using Scal = double; 
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

class KernelFactory {
 public:
  virtual Kernel* Make(Vars&, const MyBlockInfo&) = 0;
  virtual ~KernelFactory() {}
};





