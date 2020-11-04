// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#include <cstddef>
#include <string>

namespace sysinfo {

// Returns resident memory [bytes]
size_t GetMem();
// Returns true if Intel Hyperthreads are enabled
bool HasHyperthreads();

std::string GetHostname();

} // namespace sysinfo
