#pragma once

#include <array>
#include <memory>

#include "Kernel.h"

std::unique_ptr<Distr> CreateLocal(
    KernelFactory& kf, int bs, Idx b, Idx p, int es, int h);
