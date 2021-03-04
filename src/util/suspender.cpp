// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "util/logger.h"

#include "suspender.h"

Suspender::Sem::Sem(Suspender& owner, const std::string& name) : owner_(owner) {
  auto& states = owner_.states_;
  auto& pos = owner_.pos_;

  if (pos == states.begin()) {
    owner_.allow_nested_ = true; // allow nested calls on first level
    owner_.depth_ = 0;
  }

  fassert(
      owner_.allow_nested_,
      owner_.GetNameSequence() +
          ": Nested calls not allowed. Use `sem.Nested()` on upper level");
  owner_.allow_nested_ = false;

  if (std::next(pos) == states.end()) {
    states.emplace_back(0, 0, name);
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
    pos->suff = suff;
    pos->suff_id = pos->target;
    ++owner_.depth_;
    return true;
  }
  return false;
}

bool Suspender::Sem::operator()(const std::string& suff) {
  owner_.allow_nested_ = false;
  return Next(suff);
}

bool Suspender::Sem::Nested(const std::string& suff) {
  owner_.allow_nested_ = true;
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

Suspender::Suspender() : allow_nested_(false) {
  states_.emplace_back(-1, -1, "");
  pos_ = states_.begin();
}

Suspender::Sem Suspender::GetSem(const std::string& name) {
  return Sem(*this, name);
}

std::string Suspender::GetNameSequence() const {
  std::string res;
  bool first = true;
  for (const auto& s : states_) {
    if (!first) {
      res += " --> ";
    } else {
      first = false;
    }
    const int cnt = s.suff_id;
    const char c0 = '0' + cnt % 10;
    const char c1 = '0' + (cnt / 10) % 10;
    res += s.name + ":" + c1 + c0 + ":" + s.suff;
  }
  return res;
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
