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
  // Prints stats to cout
  void Report();
  // Returns current dump index
  size_t GetN() { return pn_; }

 private:
  Vars& var;
  std::string pre_;
  int pn_ = -1; // dum[p] index
  double pt_ = -std::numeric_limits<double>::max(); // last dum[p] time
  double ptt_; // dum[p] target time
};

