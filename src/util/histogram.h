#pragma once

#include <mpi.h>
#include <map>
#include <stack>
#include <string>
#include <vector>

#include "metrics.h"

class Histogram {
 public:
  Histogram(
      const MPI_Comm comm, const std::string &name, const bool active = true)
      : comm_(comm), name_(name), active_(active) {}
  ~Histogram() {
    if (active_) {
      Consolidate();
    }
  }

  Histogram(const Histogram &c) = delete;
  Histogram &operator=(const Histogram &c) = delete;

  void SeedSample() {
    if (active_) {
      timers_.push(SingleTimer());
    }
  }

  void CollectSample(const std::string &name) {
    if (active_ && !timers_.empty()) {
      samples_[name].push_back(timers_.top().GetSeconds());
      timers_.pop();
    }
  }

 private:
  const MPI_Comm comm_;
  const std::string name_;
  const bool active_;
  std::map<std::string, std::vector<double>> samples_;
  std::stack<SingleTimer> timers_;

  void Consolidate();
};
