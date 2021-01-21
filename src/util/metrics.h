// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <map>
#include <stack>

#include "timer.h"

template <class Key_>
class MultiTimer {
 public:
  using Key = Key_;
  using Value = double;

  // Marks start of timer with empty key.
  void Push();
  // Marks stop of timer, appends to accumulated time for new key.
  void Pop(const Key& key);
  // Returns map of accumulated time
  const std::map<Key, Value>& GetMap() const;
  // Resets all accumulated time to zero
  void Reset();

 private:
  SingleTimer timer_;
  struct Start {
    Key key;
    Value time;
  };
  std::map<Key, Value> accum_time_;
  std::stack<Start> starts_;
};

template <class Key>
void MultiTimer<Key>::Push() {
  starts_.push({{}, timer_.GetSeconds()});
}

template <class Key>
void MultiTimer<Key>::Pop(const Key& key) {
  Start& start = starts_.top();
  start.key = key;
  accum_time_[start.key] += timer_.GetSeconds() - start.time;
  starts_.pop();
}

template <class Key>
auto MultiTimer<Key>::GetMap() const -> const std::map<Key, Value>& {
  return accum_time_;
}

template <class Key>
void MultiTimer<Key>::Reset() {
  accum_time_.clear();
}
