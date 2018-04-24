#include <iostream>
#include <cassert>

#include "dumper.h"

bool Dumper::Try(double t, double dt) {
  // requirements: 
  // * interval between dumps is at least dumpdt + dt
  // * dumpt % dtumpdt <= dt * 0.5
  auto pd = var.Double["dumpdt"]; // dum[p] [d]t

  // next target time
  ptt_ = int(std::max<double>(pt_ + pd, 0.) / pd + 0.5) * pd;
  assert(ptt_ > pt_);

  if (var.Int["output"] && 
      t >= ptt_ - dt * 0.5 &&
      pn_ < var.Int["dumpmax"]) {
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
