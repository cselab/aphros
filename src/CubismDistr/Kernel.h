#pragma once

#include <string>

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


struct MyBlockInfo {
  int index[3];
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


using Idx = std::array<int, 3>;

template <int n=3>
Idx GetIdx(const int* d) {
  return {d[0], d[1], d[2]};
}


class Distr {
 public:
  virtual bool IsDone() const = 0;
  virtual void Step() = 0;
  virtual ~Distr() {}
};

