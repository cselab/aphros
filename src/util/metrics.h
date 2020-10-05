// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <chrono>
#include <map>
#include <stack>

template <class Key_>
class MultiTimer {
 public:
  using Key = Key_;
  using Value = double;

  // Marks start of timer with given key
  void Push(const Key& key);
  void Push() {
    Push("");
  }
  // Marks stop of timer with last key.
  // Overwrites the key if key != "".
  // Appends to accumulated time for new key.
  void Pop(const Key& key);
  void Pop() {
    Pop("");
  }
  // Returns map of accumulated time
  const std::map<Key, Value>& GetMap() const;
  // Resets all accumulated time to zero
  void Reset();

 private:
  using Clock = std::chrono::steady_clock;
  Clock clock_;
  struct Start { // start
    Key key;
    Clock::time_point time;
    Start(const Key& key_, Clock::time_point time_) : key(key_), time(time_) {}
  };
  std::map<Key, Value> accum_; // accumulated time
  std::stack<Start> starts_;
};

template <class Key>
void MultiTimer<Key>::Push(const Key& key) {
  starts_.emplace(key, clock_.now());
}

template <class Key>
void MultiTimer<Key>::Pop(const Key& key) {
  Start& start = starts_.top();
  if (key != "") {
    start.key = key;
  }
  if (!accum_.count(start.key)) {
    accum_[start.key] = 0.;
  }
  accum_[start.key] +=
      std::chrono::duration<Value>(clock_.now() - start.time).count();
  starts_.pop();
}

template <class Key>
auto MultiTimer<Key>::GetMap() const -> const std::map<Key, Value>& {
  return accum_;
}

template <class Key>
void MultiTimer<Key>::Reset() {
  accum_.clear();
}
