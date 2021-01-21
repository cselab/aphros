// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>
#include <string>
#include <utility>

class ExecutionTimer {
 public:
  // name: returned by GetName()
  // timeout: time in seconds for which to repeat F()
  // batch: number of calls F() at every iteration
  ExecutionTimer(std::string name, double timeout = 0.01, size_t batch = 1);
  ExecutionTimer(std::string name);
  virtual ~ExecutionTimer() = default;
  std::string GetName() const;
  // Repeats batches of  F() until reaching the timeout.
  // Returns minimal execution time per call and number of calls.
  struct Result {
    double min_call_time;
    size_t iters;
  };
  Result Run();
  // Function to evaluate.
  // Implementation note: use a volatile variable to prevent optimization.
  virtual void F() = 0;

 private:
  void Batch();

  const std::string name_;
  const double timeout_;
  const size_t batch_;
};

class SingleTimer {
 public:
  SingleTimer();
  SingleTimer(const SingleTimer&) = delete;
  SingleTimer(SingleTimer&&) = delete;
  SingleTimer& operator=(const SingleTimer&) = delete;
  SingleTimer& operator=(SingleTimer&&) = delete;
  ~SingleTimer();
  double GetSeconds() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
