// Created by Petr Karnakov on 11.05.2022
// Copyright 2022 ETH Zurich

#include <stdint.h>
#include <stdlib.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>

#include "util/logger.h"
#include "parse/parser.h"

Vars g_var;
std::string g_var_str;
Parser g_parser(g_var);

static void main_loop() {
  EM_ASM_({ Draw(); });
}

extern "C" {
int SetConfig(const char* str) {
  std::stringstream buf(str);
  Vars var;
  try {
    Parser(var).ParseStream(buf);
  }
  catch(const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  g_var = var;
  return 0;
}
const char* GetConfig() {
  std::stringstream buf;
  Parser(g_var).PrintVars(buf);
  g_var_str = buf.str();
  return g_var_str.c_str();
}
}

int main() {
  emscripten_set_main_loop(main_loop, 1, 0);
}
