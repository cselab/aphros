// Created by Petr Karnakov on 23.07.2020
// Copyright 2020 ETH Zurich

#include <string>

namespace util {

// Returns resolved path.
std::string GetRealpath(std::string path);

// Creates directory.
// path: path to target directory
// parent: if true, make parent directories as needed
void Makedir(std::string path, bool parent=true);

// Returns true if path is a directory.
bool IsDir(std::string path);

} // namespace util
