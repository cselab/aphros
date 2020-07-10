// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include <string>
#include <stdexcept>
#include <sstream>

#define FILELINE \
  (std::string() + __FILE__ + ":" + std::to_string(__LINE__))

// force assert
#define fassert(x)                                                        \
  do {                                                                    \
    if (!(x)) {                                                           \
      throw std::runtime_error(FILELINE + ": assertion failed '" #x "'"); \
    }                                                                     \
  } while (0);

#define NAMEVALUE(x)                 \
  ([&]() -> std::string {            \
    std::stringstream namevalue_s;   \
    namevalue_s << #x << '=' << (x); \
    return namevalue_s.str();        \
  }())

#define fassert_equal(x, y)                                                   \
  do {                                                                        \
    if (!((x) == (y))) {                                                      \
      std::stringstream fasrteq_s;                                            \
      fasrteq_s << FILELINE << ": assertion failed, expected equal ";      \
      fasrteq_s << #x << "='" << (x) << "' and " << #y << "='" << (y) << "'"; \
      throw std::runtime_error(fasrteq_s.str());                              \
    }                                                                         \
  } while (0);
