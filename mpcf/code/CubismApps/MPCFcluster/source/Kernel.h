#pragma once

#include <string>
#include <memory>

#include "../../hydro/suspender.h"

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
  virtual std::unique_ptr<Kernel> Make(const MyBlockInfo&) = 0;
  virtual ~KernelFactory() {}
};


using Idx = std::array<int, 3>;

class Distr {
 public:
  virtual bool IsDone() const = 0;
  virtual void Step() = 0;
};

