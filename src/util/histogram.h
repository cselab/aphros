// Created by Fabian Wermelinger on 29.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <map>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

#include "timer.h"
#include "util/logger.h"
#include "util/mpi.h"

class Sampler {
 public:
  Sampler(const bool active = true) : active_(active) {}
  Sampler(const Sampler&) = delete;
  Sampler(Sampler&&) = default;
  Sampler& operator=(const Sampler&) = delete;
  Sampler& operator=(Sampler&&) = delete;
  virtual ~Sampler() {}

  using SampleMap = std::map<std::string, std::vector<double>>;

  void SeedSample() {
    if (active_) {
      timers_.emplace();
    }
  }

  void CollectSample(const std::string& name) {
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

  void Append(const Sampler& s) {
    if (active_) {
      for (const auto& x : s.samples_) {
        auto& my_s = samples_[x.first];
        const auto& their_s = x.second;
        my_s.insert(my_s.end(), their_s.begin(), their_s.end());
      }
    }
  }

  void AppendSample(const std::string& name, const double samp) {
    if (active_) {
      samples_[name].push_back(samp);
    }
  }

  void Insert(const std::string& name, const std::vector<double>& data) {
    if (active_) {
      auto it = samples_.find(name);
      if (it == samples_.end()) {
        samples_[name] = data;
      } else {
        auto& my_s = it->second;
        my_s.insert(my_s.end(), data.begin(), data.end());
      }
    }
  }

  void AddTo(const std::string& addto, const std::vector<double>& yours) {
    if (active_) {
      auto& mine = samples_.at(addto);
      fassert_equal(
          mine.size(), yours.size(), ", Add: vectors are of unequal length");
      for (size_t i = 0; i < mine.size(); ++i) {
        mine[i] += yours[i];
      }
    }
  }

  void SubtractFrom(const std::string& from, const std::vector<double>& yours) {
    if (active_) {
      auto& mine = samples_.at(from);
      fassert_equal(
          mine.size(), yours.size(), ", Add: vectors are of unequal length");
      for (size_t i = 0; i < mine.size(); ++i) {
        mine[i] -= yours[i];
      }
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
      Consolidate_();
    }
  }

  Histogram(const Histogram& c) = delete;
  Histogram& operator=(const Histogram& c) = delete;

 private:
  const MPI_Comm comm_;
  const std::string name_;

  void Consolidate_();
  void HomogenizeCollection_();
};
