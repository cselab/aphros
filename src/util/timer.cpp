// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#include <cassert>
#include <chrono>
#include <limits>

#include "timer.h"

ExecutionTimer::ExecutionTimer(std::string name, double timeout, size_t batch)
    : name_(name), timeout_(timeout), batch_(batch) {}

std::string ExecutionTimer::GetName() const {
  return name_;
}

auto ExecutionTimer::Run() -> Result {
  using Clock = std::chrono::steady_clock;
  Clock clock;
  auto startall = clock.now();
  double total_time = 0;
  double min_batch_time = std::numeric_limits<double>::max();
  size_t iter = 0;
  do {
    auto start = clock.now();
    Batch();
    auto stop = clock.now();
    iter += batch_;
    total_time = std::chrono::duration<double>(stop - startall).count();
    min_batch_time = std::min(
        min_batch_time, std::chrono::duration<double>(stop - start).count());
  } while (total_time < timeout_);

  return {min_batch_time / batch_, iter};
}

void ExecutionTimer::Batch() {
  for (size_t i = 0; i < batch_; ++i) {
    F();
  }
}

struct SingleTimer::Imp {
  using Clock = std::chrono::steady_clock;
  Clock clock_;
  Clock::time_point start_ = clock_.now();

  Imp() = default;
};

SingleTimer::SingleTimer() : imp(new Imp()) {}
SingleTimer::~SingleTimer() = default;

double SingleTimer::GetSeconds() const {
  return std::chrono::duration<double>(imp->clock_.now() - imp->start_).count();
}
