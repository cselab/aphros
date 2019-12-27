#pragma once

#include "parse/vars.h"

template <class Par>
struct ParsePar {
  Par operator()(const Vars& var);
};

