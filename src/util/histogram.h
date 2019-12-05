#pragma once

#include <mpi.h>
#include <map>
#include <stack>
#include <string>
#include <vector>

#include "metrics.h"

class Sampler {
 public:
  Sampler(const bool active = true) : active_(active) {}
  virtual ~Sampler() {}

  using SampleMap = std::map<std::string, std::vector<double>>;

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

  void PopLast(const std::string& name) {
    if (active_ && !samples_[name].empty()) {
      samples_[name].pop_back();
    }
  }

  const SampleMap& GetSamples() const {
    return samples_;
  }

 protected:
  const bool active_;
  SampleMap samples_;
  std::stack<SingleTimer> timers_;
};

class Histogram : public Sampler {
 public:
  Histogram(
      const MPI_Comm comm, const std::string& name, const bool active = true)
      : Sampler(active), comm_(comm), name_(name) {}
  ~Histogram() {
    if (active_) {
      Consolidate();
    }
  }

  Histogram(const Histogram& c) = delete;
  Histogram& operator=(const Histogram& c) = delete;

  void Append(const Sampler& s) {
    if (active_) {
      for (const auto& x : s.GetSamples()) {
        auto& my_s = samples_[x.first];
        const auto& their_s = x.second;
        my_s.insert(my_s.end(), their_s.begin(), their_s.end());
      }
    }
  }

 private:
  const MPI_Comm comm_;
  const std::string name_;

  void Consolidate();
};
