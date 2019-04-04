#pragma once

#include <string>
#include <set>

// Returns set of words from space separated string.
std::set<std::string> GetWords(std::string s) {
  std::set<std::string> r;
  std::stringstream st(s);
  st >> std::skipws;
  while (true) {
    std::string e;
    st >> e;
    if (st) {
      r.insert(e);
    } else {
      break;
    }
  }
  return r;
}

