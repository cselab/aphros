/*
 * metrics.hpp
 *
 *  Created on: Mar 2, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <map>
#include <stack>
#include <chrono>

class SingleTimer {
  using Clock =  std::chrono::steady_clock;
  Clock clock_;
  Clock::time_point start_;
 public:
  SingleTimer() :
    start_(clock_.now()) {}
  double GetSeconds() const {
    return std::chrono::duration<double>(clock_.now() - start_).count();
  }
};

template <class Attr>
class MultiTimer {
  using Clock =  std::chrono::steady_clock;
  Clock clock_;
  struct Timer {
    const Attr attr_;
    Clock::time_point start_;
    Timer(const Attr& attr, Clock::time_point start)
        : attr_(attr)
        , start_(start)
    {}
  };
  std::map<Attr, double> total_time_;
  std::stack<Timer> stack_;

 public:
  void Push(const Attr& attr) {
    stack_.push(Timer(attr, clock_.now()));
  }
  void Pop() {
    const Timer& timer = stack_.top();
    if (!total_time_.count(timer.attr_)) {
      total_time_[timer.attr_] = 0.;
    }
    total_time_[timer.attr_] +=
        std::chrono::duration<double>(clock_.now() - timer.start_).count();
    stack_.pop();
  }
  const std::map<Attr, double>& GetTotalTime() const {
    return total_time_;
  }
};
