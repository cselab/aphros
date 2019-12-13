#undef NDEBUG
#include <cassert>
#include <iostream>

#include "util/suspender.h"

#include "context.h"
#include "loop.h"
#include "other.h"
#include "simple.h"

int main() {
  simple::Simple();
  other::Test();
  loop::Test();
  context::Test();
}
