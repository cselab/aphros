#pragma once

#include <limits>
#include <string>

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

// fld: field name
// ext: extension
// t: time step
// it: iteration, skip if -1
std::string GetDumpName(std::string fld, std::string ext, size_t t, size_t it);

// fld: field name
// ext: extension
// t: time step
std::string GetDumpName(std::string fld, std::string ext, size_t t);
