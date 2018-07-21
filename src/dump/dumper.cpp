#include <iostream>
#include <cassert>
#include <string>
#include <limits>
#include <sstream>
#include <iomanip>

#include "dumper.h"

bool Dumper::Try(double t, double dt) {
  if (pt_ == t) { 
    return true;
  }
  if (double* t0 = var.Double(pre_ + "t0")) {
    if (t < *t0 - dt * 0.5) {
      return false;
    }
  }

  // requirements: 
  // * interval between dumps is at least dumpdt + dt
  // * dumpt % dtumpdt <= dt * 0.5
  auto pdt = var.Double[pre_ + "dt"]; // dum[p] [d]t

  // next target time
  if (pdt > 0.) {
    ptt_ = int(std::max<double>(pt_ + pdt, 0.) / pdt + 0.5) * pdt;
  } else {
    ptt_ = t;
  }

  if (var.Int["output"] && 
      t >= ptt_ - dt * 0.5 &&
      pn_ < var.Int[pre_ + "max"]) {
    pt_ = t;
    ++pn_;
    return true;
  }
  return false;
}

void Dumper::Report() {
  std::cerr << "Dump " 
      << "n=" << pn_
      << " t=" << pt_
      << " target=" << ptt_
      << std::endl;
}

std::string GetDumpName(std::string fld, std::string ext, int t, int it) {
  std::stringstream s;
  s << fld;
  s << "_" << std::setfill('0') << std::setw(4) << t;
  if (it != -1) {
    s << "_" << std::setfill('0') << std::setw(4) << it;
  }
  s << ext;
  return s.str();
}

std::string GetDumpName(std::string fld, std::string ext, int t) {
  return GetDumpName(fld, ext, t, -1);
}
