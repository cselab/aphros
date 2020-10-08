// Created by Petr Karnakov on 23.07.2020
// Copyright 2020 ETH Zurich

#include <limits.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "filesystem.h"
#include "logger.h"

namespace util {

std::string GetRealpath(std::string path) {
  char buf[PATH_MAX + 1];
  char* ptr = realpath(path.c_str(), buf);
  return ptr ? std::string(ptr) : "";
}

std::string GetDirname(std::string path) {
  char buf[PATH_MAX + 1];
  strcpy(buf, path.c_str());
  const char* ptr = dirname(buf);
  return std::string(ptr);
}

std::string GetBasename(std::string path) {
  char buf[PATH_MAX + 1];
  strcpy(buf, path.c_str());
  const char* ptr = basename(buf);
  return std::string(ptr);
}

std::string Join(std::string a, std::string b) {
  if (b == "") {
    return a;
  }
  if (a == "") {
    return b;
  }
  if (b[0] == '/') {
    return b;
  }
  if (a.back() == '/') {
    return a + b;
  }
  return a + '/' + b;
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

bool IsFile(std::string path) {
  struct stat info;
  return (stat(path.c_str(), &info) == 0) &&
         (info.st_mode & (S_IFREG | S_IFLNK));
}

bool IsDir(std::string path) {
  struct stat info;
  return (stat(path.c_str(), &info) == 0) && (info.st_mode & S_IFDIR);
}

} // namespace util
