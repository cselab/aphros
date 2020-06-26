// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include <string>
#include <stdexcept>

#define FILELINE \
  (std::string() + __FILE__ + ":" + std::to_string(__LINE__))

// force assert
#define fassert(x)                                                        \
  do {                                                                    \
    if (!(x)) {                                                           \
      throw std::runtime_error(FILELINE + ": assertion failed '" #x "'"); \
    }                                                                     \
  } while (0);

