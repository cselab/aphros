// Created by Sergey Litvinov on 08.03.2021
// Copyright 2021 ETH Zurich

#include "system.h"
#ifdef _WIN32
#include "system_windows.inc"
#else
#include "system_unix.inc"
#endif

int SystemIsDir(const char* path) {
  struct stat info;
  return stat(path, &info) == 0 && S_ISDIR(info.st_mode);
}

int SystemIsFile(const char* path) {
  struct stat info;
  return stat(path, &info) == 0 && (S_ISREG(info.st_mode) || S_ISLNK(info.st_mode));
}
