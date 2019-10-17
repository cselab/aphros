#undef NDEBUG
#include <cassert>
#include <iostream>

#include "util/suspender.h"

#include "simple.h"
#include "other.h"
#include "loop.h"
#include "context.h"


int main() {
  simple::Simple();
  other::Test();
  loop::Test();
  context::Test();
}
