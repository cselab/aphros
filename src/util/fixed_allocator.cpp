// Created by Petr Karnakov on 30.10.2020
// Copyright 2020 ETH Zurich

#include "fixed_allocator.h"

namespace fast_allocator_impl {

std::map<size_t, FixedAllocatorGeneric*> g_alloc_map;

// Register a series of FixedAllocator's
// with the chunk size from 8 bytes till 4 megabytes
FixedAllocatorSeries<8, 10, (1 << 20)> alloc_series_default_;

} // namespace fast_allocator_impl
