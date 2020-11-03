// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include "macros.h"

#define FILELINE (std::string() + __FILE__ + ":" + std::to_string(__LINE__))

#define NAMEVALUE(x)                 \
  ([&]() -> std::string {            \
    std::stringstream namevalue_s;   \
    namevalue_s << #x << '=' << (x); \
    return namevalue_s.str();        \
  }())

// force assert
#define fassert_1(x)                                                      \
  do {                                                                    \
    if (!(x)) {                                                           \
      throw std::runtime_error(FILELINE + ": assertion failed '" #x "'"); \
    }                                                                     \
  } while (0);

#define fassert_2(x, msg)                                      \
  do {                                                         \
    if (!(x)) {                                                \
      throw std::runtime_error(                                \
          FILELINE + ": assertion failed '" #x "'\n" + (msg)); \
    }                                                          \
  } while (0);

#define fassert(...) APHROS_XCAT(fassert##_, VA_SIZE(__VA_ARGS__))(__VA_ARGS__)

#define fassert_equal_2(x, y)                                                 \
  do {                                                                        \
    if (!((x) == (y))) {                                                      \
      std::stringstream fasrteq_s;                                            \
      fasrteq_s << FILELINE << ": assertion failed, expected equal ";         \
      fasrteq_s << #x << "='" << (x) << "' and " << #y << "='" << (y) << "'"; \
      throw std::runtime_error(fasrteq_s.str());                              \
    }                                                                         \
  } while (0);

#define fassert_equal_3(x, y, msg)                                           \
  do {                                                                       \
    if (!((x) == (y))) {                                                     \
      std::stringstream fasrteq_s;                                           \
      fasrteq_s << FILELINE << ": assertion failed, expected equal ";        \
      fasrteq_s << #x << "='" << (x) << "' and " << #y << "='" << (y) << "'" \
                << msg;                                                      \
      throw std::runtime_error(fasrteq_s.str());                             \
    }                                                                        \
  } while (0);

#define fassert_equal(...) \
  APHROS_XCAT(fassert_equal##_, VA_SIZE(__VA_ARGS__))(__VA_ARGS__)
