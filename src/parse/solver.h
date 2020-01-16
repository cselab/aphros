// Created by Petr Karnakov on 27.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include "parse/vars.h"

template <class Par>
struct ParsePar {
  Par operator()(const Vars& var);
};
