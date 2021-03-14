// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <sstream>
#include <string>

#include "macros.h"

typedef void (*aphros_ErrorHandler)(int code, const char* str);

constexpr int aphros_MaxErrorString = 65536;

extern "C" {
// Sets current error code and string.
// If the error code is non-zero, calls the error handler.
void aphros_SetError(int code, const char* str);
// Returns current error code. Zero value indicates no error.
int aphros_GetErrorCode();
// Returns current error string. Undefined in case of no error.
const char* aphros_GetErrorString();
// Sets the error handler. Null resets to the default handler.
void aphros_SetErrorHandler(aphros_ErrorHandler);
// Returns current error handler.
aphros_ErrorHandler aphros_GetErrorHandler();
// Default error handler. Performs no action.
void aphros_DefaultErrorHandler(int code, const char* str);
}

inline void aphros_SetError(int code, const std::string& str) {
  aphros_SetError(code, str.c_str());
}

#define FILELINE (std::string() + __FILE__ + ":" + std::to_string(__LINE__))

#define NAMEVALUE(x)                 \
  ([&]() -> std::string {            \
    std::stringstream namevalue_s;   \
    namevalue_s << #x << '=' << (x); \
    return namevalue_s.str();        \
  }())

// force assert
#define fassert_1(x)                                                \
  do {                                                              \
    if (!(x)) {                                                     \
      aphros_SetError(1, FILELINE + ": assertion failed '" #x "'"); \
      throw std::runtime_error(aphros_GetErrorString());            \
    }                                                               \
  } while (0);

#define fassert_2(x, msg)                                                     \
  do {                                                                        \
    if (!(x)) {                                                               \
      aphros_SetError(1, FILELINE + ": assertion failed '" #x "'\n" + (msg)); \
      throw std::runtime_error(aphros_GetErrorString());                      \
    }                                                                         \
  } while (0);

#define GET_COUNT(_1, _2, _3, N, ...) N
#define fassert(...)                 \
  APHROS_ID_1(APHROS_ID_1(GET_COUNT( \
      __VA_ARGS__, fassert_3, fassert_2, fassert_1, fassert_0))(__VA_ARGS__))
#define fassert_equal(...)                                            \
  APHROS_ID_1(APHROS_ID_1(GET_COUNT(                                  \
      __VA_ARGS__, fassert_equal_3, fassert_equal_2, fassert_equal_1, \
      fassert_equal_0))(__VA_ARGS__))

#define fassert_equal_2(x, y)                                         \
  do {                                                                \
    const auto fasrteq_x = (x);                                       \
    const auto fasrteq_y = (y);                                       \
    if (!(fasrteq_x == fasrteq_y)) {                                  \
      std::stringstream fasrteq_s;                                    \
      fasrteq_s << FILELINE << ": assertion failed, expected equal "; \
      fasrteq_s << #x << "='" << fasrteq_x << "' and " << #y << "='"  \
                << fasrteq_y << "'";                                  \
      aphros_SetError(1, fasrteq_s.str());                            \
      throw std::runtime_error(aphros_GetErrorString());              \
    }                                                                 \
  } while (0);

#define fassert_equal_3(x, y, msg)                                    \
  do {                                                                \
    const auto fasrteq_x = (x);                                       \
    const auto fasrteq_y = (y);                                       \
    if (!(fasrteq_x == fasrteq_y)) {                                  \
      std::stringstream fasrteq_s;                                    \
      fasrteq_s << FILELINE << ": assertion failed, expected equal "; \
      fasrteq_s << #x << "='" << fasrteq_x << "' and " << #y << "='"  \
                << fasrteq_y << "'" << msg;                           \
      aphros_SetError(1, fasrteq_s.str());                            \
      throw std::runtime_error(aphros_GetErrorString());              \
    }                                                                 \
  } while (0);
