// Created by Petr Karnakov on 10.10.2020
// Copyright 2020 ETH Zurich
//
#include <iomanip>

#include "logger.h"

#include "format.h"

#define fassert_match(c, seq) \
  fassert(                    \
      match(c, seq),          \
      std::string() + "Expected one of \"" + seq + "\", " + report())

namespace util {

std::string Escape(std::string s) {
  std::string res;
  for (auto c : s) {
    if (c == '\n') {
      res += "\\n";
    } else if (c == '\t') {
      res += "\\t";
    } else if (c == '\r') {
      res += "\\r";
    } else {
      res += c;
    }
  }
  return res;
}

std::string ParseFormat(std::string fmt, const std::vector<Printer>& printers) {
  enum class S {
    normal,
    curly,
    givenindex,
    colon,
    width,
    precision,
    type,
    write
  };
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
  struct Mod {
    std::string givenindex;
    std::string precision;
    std::string width;
    bool leadzero = false;
    char type = 0;
  };
  Mod mod;
  size_t autoindex = 0; // current index in strs with automatic numbering
  for (size_t i = 0; i < fmt.length(); ++i) {
    char c = fmt[i];
    auto peek = [&i, &fmt]() -> char {
      return i + 1 < fmt.length() ? fmt[i + 1] : 0;
    };
    auto report = [&]() {
      return std::string() + "got '" + c + "' at position " +
             std::to_string(i) + " in \"" + Escape(fmt) + "\"";
    };
    auto reportpeek = [&]() {
      return std::string() + "got '" + peek() + "' at position " +
             std::to_string(i) + " in \"" + Escape(fmt) + "\"";
    };
    switch (state) {
      case S::normal:
        if (c == '{') {
          if (peek() == '{') { // repeated '{{'
            res += c;
            ++i;
          } else {
            mod = Mod();
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
        fassert_match(c, "}:" + digits);
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == ':') {
          state = S::colon;
        } else { // digit
          state = S::givenindex;
          mod.givenindex += c;
        }
        break;
      case S::givenindex:
        fassert_match(c, "}:" + digits);
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == ':') {
          state = S::colon;
        } else { // digit
          mod.givenindex += c;
        }
        break;
      case S::colon:
        fassert_match(c, "}.eEfgdxX" + digits);
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == '.') {
          state = S::precision;
        } else if (match(c, "eEfgdxX")) {
          state = S::type;
          --i;
        } else { // digit
          mod.leadzero = (c == '0');
          mod.width += c;
          state = S::width;
        }
        break;
      case S::width:
        fassert_match(c, "}.eEfgdxX" + digits);
        if (c == '}') {
          state = S::write;
          --i;
        } else if (c == '.') {
          state = S::precision;
        } else if (match(c, "eEfgdxX")) {
          state = S::type;
          --i;
        } else { // digit
          mod.width += c;
        }
        break;
      case S::precision:
        fassert_match(c, digits + "eEfg");
        if (match(c, digits)) {
          mod.precision += c;
        } else {
          state = S::type;
          --i;
        }
        break;
      case S::type:
        fassert_match(c, "eEfgdxX");
        mod.type = c;
        state = S::write;
        break;
      case S::write: {
        fassert_match(c, "}");
        int index;
        if (mod.givenindex != "") {
          index = toint(mod.givenindex);
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
        if (mod.precision != "") {
          ss << std::setprecision(toint(mod.precision));
        }
        if (mod.width != "") {
          ss << std::setw(toint(mod.width));
        }
        if (mod.leadzero) {
          ss << std::setfill('0') << std::internal;
        }
        switch (mod.type) {
          case 'E':
            ss << std::uppercase << std::scientific;
            break;
          case 'e':
            ss << std::scientific;
            break;
          case 'f':
            ss << std::fixed;
            break;
          case 'X':
            ss << std::uppercase << std::hex;
            break;
          case 'x':
            ss << std::hex;
            break;
          case 'g':
          case 'd':
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
