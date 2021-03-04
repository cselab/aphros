// Created by Petr Karnakov on 25.07.2018
// Copyright 2018 ETH Zurich

#include <stdexcept>
#include <string>

#include "solver.h"
#include "util/logger.h"

std::string GetName(Step l) {
  switch (l) {
    case Step::time_curr: {
      return "time_curr";
    }
    case Step::time_prev: {
      return "time_prev";
    }
    case Step::iter_curr: {
      return "iter_curr";
    }
    case Step::iter_prev: {
      return "iter_prev";
    }
    default: {
      fassert(false, "GetName: Unknown layer");
    }
  }
}

ConvSc GetConvSc(std::string s) {
  auto l = {
      ConvSc::fou, ConvSc::cd, ConvSc::sou, ConvSc::quick, ConvSc::superbee};
  for (ConvSc sc : l) {
    if (GetName(sc) == s) {
      return sc;
    }
  }
  fassert(false, "ConvSc: invalid name=" + s);
}

std::string GetName(ConvSc sc) {
  switch (sc) {
    case ConvSc::fou:
      return "fou";
    case ConvSc::cd:
      return "cd";
    case ConvSc::sou:
      return "sou";
    case ConvSc::quick:
      return "quick";
    case ConvSc::superbee:
      return "superbee";
    default:
      fassert(false, "GetName: invalid ConvSc");
  }
  return "";
}
