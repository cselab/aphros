// Created by Petr Karnakov on 05.05.2022
// Copyright 2022 ETH Zurich

#include <geom/vect.h>
#include <util/format.h>
#include <util/logger.h>
#include <cmath>
#include <iostream>
#include <sstream>

template <class T>
T RoundDiv(T a, T b) {
  return int(std::round(a)) / int(std::round(b));
}

template <class T, size_t dim>
generic::Vect<T, dim> RoundDiv(
    generic::Vect<T, dim> a, generic::Vect<T, dim> b) {
  auto res = a;
  for (size_t i = 0; i < a.size(); ++i) {
    res[i] = RoundDiv<T>(a[i], b[i]);
  }
  return a;
}

template <class T>
T from_string(const std::string& s, int* status = nullptr) {
  T res;
  std::stringstream st(s);
  st >> res;
  if (status) {
    if (st.fail()) {
      *status = 1; // Error.
    } else if (!st.eof()) {
      *status = 2; // Trailing characters.
    } else {
      *status = 0;
    }
  }
  return res;
}

template <class T, class Iterator>
T Eval(Iterator begin, Iterator end, Iterator origin);

template <class T, class Iterator>
T EvalOp(
    Iterator lbegin, Iterator lend, char op, Iterator rbegin, Iterator rend,
    Iterator origin) {
  switch (op) {
    case '+':
      return Eval<T>(lbegin, lend, origin) + Eval<T>(rbegin, rend, origin);
    case '-':
      return Eval<T>(lbegin, lend, origin) - Eval<T>(rbegin, rend, origin);
    case '*':
      return Eval<T>(lbegin, lend, origin) * Eval<T>(rbegin, rend, origin);
    case '/':
      if (lbegin != lend && *(lend - 1) == '/') {
        return RoundDiv(
            Eval<T>(lbegin, lend - 1, origin), Eval<T>(rbegin, rend, origin));
      }
      return Eval<T>(lbegin, lend, origin) / Eval<T>(rbegin, rend, origin);
    default:
      fassert(
          false, util::Format(
                     "Unknown operation: '{}' with operands \"{}\" and \"{}\"",
                     op, std::string(lbegin, lend), std::string(rbegin, rend)));
  }
}

template <class T, class Iterator>
T Eval(Iterator begin, Iterator end, Iterator origin) {
  auto failmsg = [&]() {
    return util::Format(
        "Failed at position {:} while processing \"{}[{}]\": ", begin - origin,
        std::string(origin, begin), std::string(begin, end));
  };
  // Strip leading and trailing spaces.
  while (begin != end && *begin == ' ') {
    ++begin;
  }
  while (end != begin && *(end - 1) == ' ') {
    --end;
  }
  // Jumps back to the matching parenthesis. Returns true if `pos` is modified.
  auto skip_paren_back = [&](Iterator& pos) {
    if (*pos == ')') {
      int level = 1;
      do {
        --pos;
        fassert(
            pos != begin || *pos == '(', //
            util::Format(
                failmsg() + "Unmatched parentheses at position {:} in \"{}\"",
                pos - begin, std::string(begin, end)));
        if (*pos == ')') {
          ++level;
        } else if (*pos == '(') {
          --level;
          fassert(
              level >= 0, //
              util::Format(
                  failmsg() + "Too many opening parentheses at {:} in \"{}\"",
                  pos - begin, std::string(begin, end)));
          if (level == 0) {
            return true;
          }
        }
      } while (pos != begin);
    }
    return false;
  };
  fassert(begin != end, failmsg() + "Expression is empty");
  { // Strip parentheses enclosing the whole expression.
    Iterator pos = end - 1;
    if (skip_paren_back(pos) && pos == begin) {
      return Eval<T>(begin + 1, end - 1, origin);
    }
  }
  // Look for last addition or subtraction.
  for (Iterator pos = end; pos != begin;) {
    skip_paren_back(--pos);
    switch (*pos) {
      case '+':
      case '-':
        return EvalOp<T>(begin, pos, *pos, pos + 1, end, origin);
      default:
        break;
    }
  }
  // Otherwise, look for last multiplication or division.
  for (Iterator pos = end; pos != begin;) {
    skip_paren_back(--pos);
    switch (*pos) {
      case '*':
      case '/':
        return EvalOp<T>(begin, pos, *pos, pos + 1, end, origin);
      default:
        break;
    }
  }
  // Otherwise, convert to a number.
  int status;
  const T res = from_string<T>({begin, end}, &status);
  fassert(
      status != 2,
      util::Format(
          failmsg() + "Trailing characters while converting \"{}\" to a number",
          std::string(begin, end)));
  fassert(
      !status, util::Format(
                   failmsg() + "Cannot convert \"{}\" to a number",
                   std::string(begin, end)));
  return res;
}

template <class T, class Iterator>
T EvalExpr(Iterator begin, Iterator end) {
  return Eval<T>(begin, end, begin);
}

template <class T>
T EvalExpr(const std::string& expr) {
  return Eval<T>(expr.begin(), expr.end(), expr.begin());
}

template double EvalExpr(
    std::string::const_iterator, std::string::const_iterator);
template double EvalExpr(const char*, const char*);
template double EvalExpr(const std::string& expr);
