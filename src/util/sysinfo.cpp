// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#include <unistd.h>
#include <fstream>
#include <sstream>

#include "sysinfo.h"
#include "util/logger.h"

namespace sysinfo {

// Returns value from /proc/self/status found by key k
size_t GetProcStat(std::string k) {
  std::ifstream f("/proc/self/status");

  f >> std::skipws;
  while (f) {
    std::string s;
    f >> s;
    if (s == k + ":") {
      size_t r = 0;
      f >> r;
      return r;
    }
  }
  return 0;
}

size_t GetMem() {
  return GetProcStat("VmRSS") << 10;
}

bool HasHyperthreads() {
  bool has_ht = false;
  std::ifstream f("/proc/cpuinfo");
  f >> std::skipws;
  while (f) {
    std::string s;
    f >> s;
    if (s == "flags") {
      break;
    }
  }
  std::string line;
  std::getline(f, line);
  std::istringstream is(line);
  while (!is.eof()) {
    std::string s;
    is >> s;
    if (s == "ht") {
      has_ht = true;
      break;
    }
  }
  return has_ht;
}

std::string GetHostname() {
  const size_t kMaxLength = 4096;
  char buf[kMaxLength];
  int err = gethostname(buf, kMaxLength);
  fassert_equal(err, 0);
  return buf;
}

} // namespace sysinfo
