// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include "dumper.h"

Dumper::Dumper(const Vars& var_, std::string pre) : var(var_), prefix_(pre) {}

bool Dumper::Try(double t, double dt) {
  if (time_ == t) {
    return true;
  }
  if (const double* t0 = var.Double.Find(prefix_ + "t0")) {
    if (t < *t0 - dt * 0.5) {
      return false;
    }
  }

  // requirements:
  // * interval between dumps is at least dumpdt + dt
  // * dumpt % dtumpdt <= dt * 0.5
  auto dumpdt = var.Double[prefix_ + "dt"]; // dum[p] [d]t

  // next target time
  if (dumpdt > 0.) {
    target_ = int(std::max<double>(time_ + dumpdt, 0.) / dumpdt + 0.5) * dumpdt;
  } else {
    target_ = t;
  }

  if (var.Int["output"] && t >= target_ - dt * 0.5 &&
      frame_ < var.Int[prefix_ + "max"]) {
    time_ = t;
    ++frame_;
    return true;
  }
  return false;
}

void Dumper::Report(std::ostream& out) {
  out << "Dump "
      << "n=" << frame_ << " t=" << time_ << " target=" << target_ << std::endl;
}

std::string GetDumpName(std::string fld, std::string ext, int t, int iter) {
  std::stringstream s;
  s << fld;
  s << "_" << std::setfill('0') << std::setw(4) << t;
  if (iter != -1) {
    s << "_" << std::setfill('0') << std::setw(4) << iter;
  }
  s << ext;
  return s.str();
}

std::string GetDumpName(std::string fld, std::string ext, int t) {
  return GetDumpName(fld, ext, t, -1);
}
