#include <string>
#include <fstream>
#include <sstream>

#include "sysinfo.h"

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

} // namespace sysinfo
