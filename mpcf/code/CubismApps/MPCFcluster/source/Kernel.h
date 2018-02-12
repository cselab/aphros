#pragma once

#include <string>
#include <memory>

#include "../../hydro/suspender.h"
//#include "BlockInfo.h"

struct MyBlockInfo {
  int index[3];
  void* ptrBlock;
  double h_gridpoint[3];
  double origin[3];
};

// Suspendable kernel
class Kernel {
 public:
  virtual void Run() = 0;
  //virtual void ReadBuffer(LabMPI&) = 0;
  //virtual void WriteBuffer(Block_t&) = 0;
  virtual ~Kernel() {}
  bool Pending() const {
    return susp_.Pending();
  }
  using Sem = Suspender::Sem;
  // Create semaphore (see Suspender)
  Sem GetSem(std::string name="") {
    return susp_.GetSem(name);
  }
  
 private:
  Suspender susp_;
};

class KernelFactory {
 public:
  virtual std::unique_ptr<Kernel> Make(const MyBlockInfo&) = 0;
  virtual ~KernelFactory() {}
};

