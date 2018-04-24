#pragma once

#include "parse/vars.h"

// TODO: rename to e.g. DumpScheduler


class Dumper {
 public:
  Dumper(Vars& var) : var(var) {}
  bool Try(double t,  // current time
           double dt // simulation time step
           );
  // Print stats to cout
  void Report();

 private:
  Vars& var;
  int pn_ = -1; // dum[p] index
  double pt_ = -1e10; // dum[p] time
  double ptt_; // dum[p] target time
};

