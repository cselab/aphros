// Created by Petr Karnakov on 04.03.2021
// Copyright 2021 ETH Zurich

#include <cstdio>
#include <cstring>
#include <stdexcept>

#include "logger.h"

static int error_code;
static char error_str[aphros_MaxErrorString + 1];
aphros_ErrorHandler error_handler = aphros_DefaultErrorHandler;

void aphros_SetError(int code, const char* str) {
  error_code = code;
  if (code) {
    std::strncpy(error_str, str, aphros_MaxErrorString);
    error_str[aphros_MaxErrorString] = '\0';
  }
  error_handler(code, str);
}
int aphros_GetErrorCode() {
  return error_code;
}
const char* aphros_GetErrorString() {
  return error_str;
}
void aphros_SetErrorHandler(aphros_ErrorHandler handler) {
  if (handler) {
    error_handler = handler;
  } else {
    error_handler = aphros_DefaultErrorHandler;
  }
}
aphros_ErrorHandler aphros_GetErrorHandler() {
  return error_handler;
}
void aphros_DefaultErrorHandler(int, const char*) {}
