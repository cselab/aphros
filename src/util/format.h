// Created by Petr Karnakov on 10.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace util {

inline void AppendStrings(std::vector<std::string>&) {}

// Converts `args` to strings and appends to `strs`.
template <class T, class... Args>
void AppendStrings(
    std::vector<std::string>& strs, const T& value, const Args&... args) {
  std::stringstream s;
  s << value;
  strs.push_back(s.str());
  AppendStrings(strs, args...);
}

template <class... Args>
std::vector<std::string> GetStrings(const Args&... args) {
  std::vector<std::string> res;
  AppendStrings(res, args...);
  return res;
}

using Printer = std::function<void(std::ostream&)>;

inline void AppendPrinters(std::vector<Printer>&) {}

template <class T, class... Args>
void AppendPrinters(
    std::vector<Printer>& printers, const T& value, const Args&... args) {
  printers.push_back([value](std::ostream& out) { out << value; });
  AppendPrinters(printers, args...);
}

template <class... Args>
std::vector<Printer> GetPrinters(const Args&... args) {
  std::vector<Printer> res;
  AppendPrinters(res, args...);
  return res;
}

// Parses format string and replaces numbered substitutions with given strings.
// fmt: format string
// printers: printers from arguments
// Example:
//   ParseFormat("{0} {1} {} {2} {}", {"a", "b", "c"}) == "a b a c b"
std::string ParseFormat(std::string fmt, const std::vector<Printer>& printers);

template <class... Args>
std::string Format(std::string fmt, const Args&... args) {
  auto printers = GetPrinters(args...);
  return ParseFormat(fmt, printers);
}

} // namespace util
