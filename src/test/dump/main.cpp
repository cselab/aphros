// Created by Petr Karnakov on 17.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>

#include "dump/dumper.h"

void TestName() {
  assert(GetDumpName("vf", ".h5", 8, 100) == "vf_0008_0100.h5");
  assert(GetDumpName("vf", ".h5", 8) == "vf_0008.h5");
}

/*
void TestHdf() {
  const std::string hdfpath = "test.h5";
  const std::string dname = "q";
  if (sem.Nested()) {
    Hdf<M>::Read(
        const_cast<FieldCell<Scal>&>(as_->GetField()), "testread.h5", m, dname);
  }
  if (sem.Nested()) {
    Hdf<M>::Write(as_->GetField(), hdfpath, m, dname);
  }
  if (sem()) {
    Hdf<M>::WriteXmf("test.xmf", "test", hdfpath, m, dname);
  }
}
*/

int main() {
  TestName();
}
