// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#include <cassert>
#include <chrono>

#include "timer.h"

Timer::Timer(std::string name, double timeout /*sec*/, size_t batch)
    : n_(name), to_(timeout), b_(batch) {}

Timer::Timer(std::string name, double timeout) : Timer(name, timeout, 1) {}

Timer::Timer(std::string name) : Timer(name, 0.01, 1) {}

std::string Timer::GetName() const {
  return n_;
}

std::pair<double, size_t> Timer::Run() {
  using C = std::chrono::steady_clock; // clock
  C c;
  C::time_point s = c.now(); // start
  double e; // execution time in sec
  size_t i = 0;
  do {
    B();
    i += b_;
    e = std::chrono::duration<double>(c.now() - s).count();
  } while (e < to_);

  assert(i > 0);

  return std::make_pair(e / i, i);
}

void Timer::B() {
  for (size_t i = 0; i < b_; ++i) {
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
