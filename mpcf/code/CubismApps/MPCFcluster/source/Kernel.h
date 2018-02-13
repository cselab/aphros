#pragma once

#include <string>
#include <memory>

#include "../../hydro/suspender.h"
//#include "BlockInfo.h"

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

