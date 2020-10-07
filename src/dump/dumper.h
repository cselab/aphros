// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <iosfwd>
#include <limits>
#include <string>

#include "parse/vars.h"

// TODO: rename to e.g. DumpScheduler

class Dumper {
 public:
  // var: config
  // pre: prefix for variables
  Dumper(const Vars& var, std::string pre);
  // t: current time
  // dt: simulation time step (to select dump time closest to the target)
  bool Try(double t, double dt);
  // Prints stats to cout
  void Report(std::ostream& out);
  // Returns current dump index
  size_t GetN() {
    return frame_;
  }

 private:
  const Vars& var;
  std::string prefix_;
  int frame_ = -1;
  double time_ = -std::numeric_limits<double>::max(); // last dump time
  double target_ = 0; // last dump target time
};

// fld: field name
// ext: extension
// t: time step
// iter: iteration, ignore if -1
std::string GetDumpName(std::string fld, std::string ext, int t, int iter);

// fld: field name
// ext: extension
// t: time step
std::string GetDumpName(std::string fld, std::string ext, int t);
