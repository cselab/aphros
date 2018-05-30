#include <iostream>
#include <cassert>
#include <string>
#include <limits>

#include "dumper.h"

bool Dumper::Try(double t, double dt) {
  if (pt_ == t) { 
    return true;
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
