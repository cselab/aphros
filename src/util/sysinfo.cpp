#include <string>
#include <fstream>

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

} // namespace sysinfo 
