// Created by Petr Karnakov on 10.10.2020
// Copyright 2020 ETH Zurich

#include "logger.h"

#include "format.h"

namespace util {

std::string ParseFormat(
    std::string fmt, const std::vector<std::string>& strs) {
  enum class S { normal, curly };
  auto toint = [&fmt](std::string str) {
    std::stringstream ss(str);
    int a = 0;
    ss >> a;
    fassert(
        !ss.fail(),
        FILELINE + ": Can't parse '" + str + "' as int in '" + fmt + "'");
    return a;
  };
  S state = S::normal;
  std::string res;
  std::string buf;
  size_t kauto = 0; // current index in strs with automatic numbering
  for (size_t i = 0; i < fmt.length(); ++i) {
    char c = fmt[i];
    auto peek = [&i, &fmt]() -> char {
      return i + 1 < fmt.length() ? fmt[i + 1] : 0;
    };
    switch (state) {
      case S::normal:
        if (c == '{') {
          if (peek() == '{') {
            res += c;
            ++i;
          } else {
            buf = "";
            state = S::curly;
          }
        } else if (c == '}') {
          fassert(
              peek() == '}', //
              FILELINE + ": Expected double }}, got '" + peek() + "' in '" +
                  fmt + "' at position " + std::to_string(i));
          res += c;
          ++i;
        } else {
          res += c;
        }
        break;
      case S::curly:
        if (c == '}') {
          if (buf == "") {
            fassert(
                kauto < strs.size(), //
                FILELINE + ": More fields {} in '" + fmt + "' than " +
                    std::to_string(strs.size()) + " arguments");
            res += strs[kauto];
            ++kauto;
          } else {
            auto k = toint(buf);
            fassert(
                k >= 0 && k < (int)strs.size(), //
                FILELINE + ": Invalid numbered field {" + std::to_string(k) +
                    "} in '" + fmt + "' for " + std::to_string(strs.size()) +
                    " provided arguments");
            res += strs[k];
          }
          state = S::normal;
        } else {
          buf += c;
        }
        break;
    }
  }
  return res;
}

} // namespace util
