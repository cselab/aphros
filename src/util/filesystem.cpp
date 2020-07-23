// Created by Petr Karnakov on 23.07.2020
// Copyright 2020 ETH Zurich

#include <cstdlib>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdexcept>

#include "logger.h"
#include "filesystem.h"

namespace util {

std::string GetRealpath(std::string path) {
  char buf[PATH_MAX + 1];
  char* ptr = realpath(path.c_str(), buf);
  return std::string(ptr);
}

void Makedir(std::string path, bool parent) {
  std::string cmd;
  cmd += "mkdir ";
  if (parent) {
    cmd += "-p ";
  }
  cmd += "'" + path + "'";
  std::system(cmd.c_str());

  struct stat info;
  fassert(stat(path.c_str(), &info) == 0, "Can't access '" + path + "'");
  fassert(info.st_mode & S_IFDIR, "Not a directory '" + path + "'");
}

bool IsDir(std::string path) {
  struct stat info;
  return (stat(path.c_str(), &info) == 0) && (info.st_mode & S_IFDIR);
}

} // namespace util
