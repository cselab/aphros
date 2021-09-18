// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// Returns set of words from space separated string.
inline std::unordered_set<std::string> GetWords(std::string s) {
  std::unordered_set<std::string> r;
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

// Splits string by delimiter.
inline std::vector<std::string> SplitByDelimiter(std::string str, char delim) {
  std::vector<std::string> ss;
  std::istringstream f(str);
  std::string s;
  while (std::getline(f, s, delim)) {
    ss.push_back(s);
  }
  return ss;
}
