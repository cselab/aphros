// Created by Petr Karnakov on 23.07.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <string>

namespace util {

// Returns resolved path.
std::string GetRealpath(std::string path);

std::string GetDirname(std::string path);

std::string GetBasename(std::string path);

std::array<std::string, 2> SplitExt(std::string path);

std::string Join(std::string path0, std::string path1);

// Creates directory.
// path: path to target directory
// parent: if true, make parent directories as needed
void Makedir(std::string path, bool parent = true);

// Returns true if path is a file.
bool IsFile(std::string path);

// Returns true if path is a directory.
bool IsDir(std::string path);

} // namespace util
