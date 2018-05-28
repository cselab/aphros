#pragma once

#include "parse/vars.h"

// TODO: rename to e.g. DumpScheduler


class Dumper {
 public:
  // var: config
  // pre: prefix for variables
  Dumper(Vars& var, std::string pre) : var(var), pre_(pre) {}
  bool Try(double t /*current time*/, 
           double dt /*simulation time step*/);
  // Print stats to cout
  void Report();

 private:
  Vars& var;
  std::string pre_;
  int pn_ = -1; // dum[p] index
  double pt_ = std::numeric_limits<double>::min(); // last dum[p] time
  double ptt_; // dum[p] target time
};

