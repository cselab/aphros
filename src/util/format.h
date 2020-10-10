#include <iostream>
#include <sstream>
#include <vector>

#include "logger.h"

namespace util {

void AppendStrings(std::vector<std::string>&) {}

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
std::string ParseFormat(std::string fmt, const std::vector<std::string>& strs) {
  enum class S { normal, curly };
  auto toint = [](std::string str) {
    std::stringstream ss(str);
    int a = 0;
    ss >> a;
    fassert(!ss.fail(), FILELINE + "Can't parse '" + str + "' as int");
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
              FILELINE + "Expected double }}, got '" + peek() + "' in '" + fmt +
                  "' at position " + std::to_string(i));
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
                FILELINE + "More fields {} than " +
                    std::to_string(strs.size()) + " arguments");
            res += strs[kauto];
            ++kauto;
          } else {
            auto k = toint(buf);
            fassert(
                k >= 0 && k < (int)strs.size(), //
                FILELINE + "Invalid numbered field {" + std::to_string(k) +
                    "} for " + std::to_string(strs.size()) +
                    " provided arguments.");
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

template <class... Args>
std::string Format(std::string fmt, Args... args) {
  std::vector<std::string> strs;
  AppendStrings(strs, args...);
  return ParseFormat(fmt, strs);
}

} // namespace util
