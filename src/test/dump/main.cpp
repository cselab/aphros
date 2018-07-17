#undef NDEBUG
#include <cassert>
#include <iostream>

#include "dump/dumper.h"

void TestName() {
  assert(GetDumpName("vf", ".h5", 8, 100) == "vf_0008_0100.h5");
  assert(GetDumpName("vf", ".h5", 8) == "vf_0008.h5");
}

int main() {
  TestName();
}
