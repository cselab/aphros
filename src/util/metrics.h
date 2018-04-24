#pragma once

#include <map>
#include <stack>
#include <chrono>

class SingleTimer {
  using Clock =  std::chrono::steady_clock;
  Clock clock_;
  Clock::time_point start_;
 public:
  SingleTimer() : start_(clock_.now()) {}
  double GetSeconds() const {
    return std::chrono::duration<double>(clock_.now() - start_).count();
  }
};

template <class Key_>
class MultiTimer {
 public:
  using Key = Key_;
  using Value = double;

  // Marks start of timer with given key
  void Push(const Key& k);
  void Push() { Push(""); }
  // Marks stop of timer with last key.
  // Overwrites the key if k != "".
  // Appends to accumulated time for new key.
  void Pop(const Key& k);
  void Pop() { Pop(""); }
  // Returns map of accumulated time
  const std::map<Key, Value>& GetMap() const;

 private:
  using Clock = std::chrono::steady_clock;
  Clock c_;
  struct S { // start 
    Key k;
    Clock::time_point t;
    S(const Key& k, Clock::time_point t) : k(k), t(t) {}
  };
  std::map<Key, Value> a_; // accumulated time
  std::stack<S> ss_;        // stack of starts
};

template <class K>
void MultiTimer<K>::Push(const Key& k) {
  ss_.emplace(k, c_.now());
}

template <class K>
void MultiTimer<K>::Pop(const Key& k) {
  S& s = ss_.top();
  if (k != "") {
    s.k = k;   
  }
  if (!a_.count(s.k)) {
    a_[s.k] = 0.;
  }
  a_[s.k] += std::chrono::duration<Value>(c_.now() - s.t).count();
  ss_.pop();
}

template <class K>
auto MultiTimer<K>::GetMap() const -> const std::map<Key, Value>&  {
  return a_;
}
