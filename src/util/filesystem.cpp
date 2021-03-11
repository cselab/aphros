// Created by Petr Karnakov on 23.07.2020
// Copyright 2020 ETH Zurich

#include <limits.h>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "system.h"
#include "filesystem.h"
#include "logger.h"

namespace util {

std::string GetRealpath(std::string path) {
  char buf[PATH_MAX + 1];
  SystemRealPath(path.c_str(), buf);
  return buf[0] != '\0' ? std::string(buf) : "";
}

std::string GetDirname(std::string path) {
  char buf[PATH_MAX + 1];
  SystemDirName(path.c_str(), buf);
  return std::string(buf);
}

std::string GetBasename(std::string path) {
  char buf[PATH_MAX + 1];
  SystemBaseName(path.c_str(), buf);
  return std::string(buf);
}

std::array<std::string, 2> SplitExt(std::string path) {
  const size_t i = path.find_last_of("/.");
  if (i == std::string::npos || path[i] == '/' || // no period in filename
      i == 0 || path[i - 1] == '/' // filename starts with period
  ) {
    return {path, ""};
  }
  return {path.substr(0, i), path.substr(i, std::string::npos)};
}

std::string Join(std::string a, std::string b) {
  char buf[PATH_MAX + 1];
  SystemJoin(a.c_str(), b.c_str(), buf);
  return buf;
}

void Makedir(std::string path, bool parent) {
  if (SystemMakeDir(path.c_str(), parent) != 0)
    fassert(0, "Can't access '" + path + "'");
}

bool IsFile(std::string path) {
  return SystemIsFile(path.c_str());
}

bool IsDir(std::string path) {
  return SystemIsDir(path.c_str());
}

} // namespace util
