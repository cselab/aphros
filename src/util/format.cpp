// Created by Petr Karnakov on 10.10.2020
// Copyright 2020 ETH Zurich
//
#include <iomanip>

#include "logger.h"

#include "format.h"

namespace util {

std::string ParseFormat(std::string fmt, const std::vector<Printer>& printers) {
  enum class S { normal, curly, givenindex, colon, precision, type, write };
  auto toint = [&fmt](std::string str) {
    std::stringstream ss(str);
    int a = 0;
    ss >> a;
    fassert(
        !ss.fail(),
        FILELINE + ": Can't parse '" + str + "' as int in '" + fmt + "'");
    return a;
  };
  auto match = [](char c, std::string s) {
    return s.find(c) != std::string::npos;
  };
  const std::string digits = "0123456789";
  S state = S::normal;
  std::string res;
  std::string givenindex;
  std::string precision;
  char type;
  size_t autoindex = 0; // current index in strs with automatic numbering
  for (size_t i = 0; i < fmt.length(); ++i) {
    char c = fmt[i];
    auto peek = [&i, &fmt]() -> char {
      return i + 1 < fmt.length() ? fmt[i + 1] : 0;
    };
    auto report = [&]() {
      return std::string() + "got '" + c + "' at position " +
             std::to_string(i) + " in '" + fmt + "'";
    };
    auto reportpeek = [&]() {
      return std::string() + "got '" + peek() + "' at position " +
             std::to_string(i) + " in '" + fmt + "'";
    };
    switch (state) {
      case S::normal:
        if (c == '{') {
          if (peek() == '{') { // repeated '{{'
            res += c;
            ++i;
          } else {
            givenindex = "";
            precision = "";
            type = 0;
            state = S::curly;
          }
        } else if (c == '}') {
          fassert(
              peek() == '}',
              FILELINE + ": Expected double }}, " + reportpeek());
          res += c;
          ++i;
        } else {
          res += c;
        }
        break;
      case S::curly:
        fassert(match(c, "}:" + digits), report());
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == ':') {
          state = S::colon;
        } else { // digit
          state = S::givenindex;
          givenindex += c;
        }
        break;
      case S::givenindex:
        fassert(match(c, "}:" + digits), report());
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == ':') {
          state = S::colon;
        } else { // digit
          givenindex += c;
        }
        break;
      case S::colon:
        fassert(match(c, "}.efg"), report());
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == '.') {
          state = S::precision;
        } else {
          state = S::type;
          --i;
        }
        break;
      case S::precision:
        fassert(match(c, digits + "efg"), report());
        if (match(c, digits)) {
          precision += c;
        } else {
          state = S::type;
          --i;
        }
        break;
      case S::type:
        fassert(match(c, "efg"), report());
        type = c;
        state = S::write;
        break;
      case S::write: {
        fassert(match(c, "}"), report());
        int index;
        if (givenindex != "") {
          index = toint(givenindex);
          fassert(
              index >= 0 && index < (int)printers.size(), //
              FILELINE + ": Invalid numbered field {" + std::to_string(index) +
                  "} in '" + fmt + "' for " + std::to_string(printers.size()) +
                  " provided arguments; " + report());
        } else {
          fassert(
              autoindex < printers.size(), //
              FILELINE + ": More fields {} in '" + fmt + "' than " +
                  std::to_string(printers.size()) + " arguments; " + report());
          index = autoindex;
          ++autoindex;
        }
        std::stringstream ss;
        if (precision != "") {
          ss << std::setprecision(toint(precision));
        }
        switch (type) {
          case 'e':
            ss << std::scientific;
            break;
          case 'f':
            ss << std::fixed;
            break;
          case 'g':
          default:
            break;
        }
        printers[index](ss);
        res += ss.str();
        state = S::normal;
        break;
      }
    }
  }
  return res;
}

} // namespace util
