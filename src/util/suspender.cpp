// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "suspender.h"

Suspender::Sem::Sem(Suspender& owner, const std::string& name)
    : owner_(owner), name_(name) {
  auto& states = owner_.states_;
  auto& pos = owner_.pos_;

  if (pos == states.begin()) {
    owner_.nest_ = true; // allow nested calls on first level
    // owner_.curname_ = "";
    owner_.depth_ = 0;
  }

  if (!owner_.nest_) {
    throw std::runtime_error(
        owner_.GetCurName() +
        ": Nested calls not allowed. Use Nested() on upper level");
  }
  owner_.nest_ = false;

  if (std::next(pos) == states.end()) {
    states.emplace_back(0, 0);
  }
  ++pos;

  pos->current = 0;
}

Suspender::Sem::~Sem() {
  auto& states = owner_.states_;
  auto& pos = owner_.pos_;

  assert(!states.empty());
  assert(pos != states.end());
  assert(pos != states.begin());

  auto ip = std::prev(pos);

  if (std::next(pos) == states.end()) {
    // all lower levels done, next stage
    ++pos->target;
    // pos->current keeps total number of stages
    if (pos->current == pos->target || pos->current == 0) {
      // all stages done or no stages, remove current level
      states.pop_back();
    }
  }
  pos = ip;
}

bool Suspender::Sem::Next(const std::string& suff) {
  auto& pos = owner_.pos_;
  if (pos->current++ == pos->target) {
    /*
    if (owner_.curname_ != "") {
      owner_.curname_ += " --> ";
    }
    std::stringstream st;
    st << std::setfill('0') << std::setw(2) << pos->target;
    owner_.curname_ += name_ + ":" + st.str() + ":" + suff;
    */
    ++owner_.depth_;
    return true;
  }
  return false;
}

bool Suspender::Sem::operator()(const std::string& suff) {
  owner_.nest_ = false;
  return Next(suff);
}

bool Suspender::Sem::Nested(const std::string& suff) {
  owner_.nest_ = true;
  return Next(suff);
}

void Suspender::Sem::LoopBegin() {
  State& s = *owner_.pos_;
  if (s.current == s.target) {
    if (s.loop_begin < s.current) {
      s.loop_begin = s.current;
    }
    ++s.target;
  }
  ++s.current;
}

void Suspender::Sem::LoopBreak() {
  State& s = *owner_.pos_;
  s.target = s.loop_end; // set target beyond loop end
}

// Important to initialize loop_end even before target reaches LoopEnd
// to be able to break on first iteration
void Suspender::Sem::LoopEnd() {
  State& s = *owner_.pos_;
  if (s.loop_end < s.loop_begin && // loop_end not initialized for current loop
      s.loop_begin <= s.target) { // target is within the loop
    s.loop_end = s.current;
  }
  if (s.current == s.target + 1) { // passed loop end
    s.target = s.loop_begin; // reset target to loop begin
  }
  ++s.current;
}

Suspender::Suspender() : nest_(false) {
  states_.emplace_back(-1, -1);
  pos_ = states_.begin();
}

Suspender::Sem Suspender::GetSem(const std::string& name) {
  return Sem(*this, name);
}

const std::string& Suspender::GetCurName() const {
  return curname_;
}

std::string Suspender::Print() const {
  std::stringstream b;
  for (auto& e : states_) {
    b << "(" << e.current << " " << e.target << ") ";
  }
  return b.str();
}

bool Suspender::Pending() const {
  return states_.size() != 1;
}
