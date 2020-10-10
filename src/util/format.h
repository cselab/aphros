// Created by Petr Karnakov on 10.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <sstream>
#include <string>
#include <vector>

namespace util {

inline void AppendStrings(std::vector<std::string>&) {}

// Converts `args` to strings and appends to `strs`.
template <class T, class... Args>
void AppendStrings(
    std::vector<std::string>& strs, const T& value, Args... args) {
  std::stringstream s;
  s << value;
  strs.push_back(s.str());
  AppendStrings(strs, args...);
}

// Parses format string and replaces numbered substitutions with given strings.
// fmt: format string
// strs: strings to substitute with
// Example:
//   ParseFormat("{0} {1} {} {2} {}", {"a", "b", "c"}) == "a b a c b"
std::string ParseFormat(std::string fmt, const std::vector<std::string>& strs);

template <class... Args>
std::string Format(std::string fmt, Args... args) {
  std::vector<std::string> strs;
  AppendStrings(strs, args...);
  return ParseFormat(fmt, strs);
}

} // namespace util
