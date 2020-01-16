// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <set>
#include <string>
#include <vector>

// Returns set of words from space separated string.
inline std::set<std::string> GetWords(std::string s) {
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

// Splits string by separator
inline std::vector<std::string> Split(std::string str, char sep) {
  std::vector<std::string> ss;
  std::istringstream f(str);
  std::string s;
  while (std::getline(f, s, sep)) {
    ss.push_back(s);
  }
  return ss;
}
