// Created by Petr Karnakov on 16.02.2018
// Copyright 2018 ETH Zurich

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
