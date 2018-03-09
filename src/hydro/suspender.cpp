#include <cassert>

#include "suspender.h"

Suspender::Sem::Sem(Suspender& p, std::string name) 
  : p(p), name_(name)
{
  auto& l = p.lu_;
  auto& i = p.lui_;

  if (i == l.begin()) {
    // Allow nested calls on first level
    p.nest_ = true;
    // Clear curname_
    p.curname_ = "";
    p.depth_ = 0;
  }

  assert(p.nest_ && "Nested calls not allowed. Use Nested() on upper level");
  p.nest_ = false;

  if (std::next(i) == l.end()) {
    l.emplace_back(0, 0);
  }
  ++i;

  i->c = 0;

}

Suspender::Sem::~Sem() {
  auto& l = p.lu_;
  auto& i = p.lui_;

  assert(!l.empty());
  assert(i != l.end());
  assert(i != l.begin());

  auto ip = std::prev(i);

  if (std::next(i) == l.end()) {
    // all lower levels done, next stage
    ++i->t;
    // i->c keeps total number of stages
    if (i->c == i->t || i->c == 0) { 
      // all stages done or no stages, remove current level
      l.pop_back();
    }
  } 
  i = ip;
}

bool Suspender::Sem::operator()(std::string suff) {
  auto& i = p.lui_;
  if (i->c++ == i->t) {
    if (p.curname_ != "") {
      p.curname_ += " --> ";
    }
    p.curname_ += name_ + ":" + suff;
    ++p.depth_;
    return true;
  }
  return false;
}

bool Suspender::Sem::Nested(std::string suff) {
  p.nest_ = true;
  return (*this)(suff);
}

#include <iostream>
void Suspender::Sem::LoopBegin() {
  if ((*this)()) {
    U& u = *p.lui_;
    u.lb = u.c - 1;
    std::cerr << "lb=" << u.lb << std::endl;
  }
}

void Suspender::Sem::LoopBreak() {
  if ((*this)()) {
    U& u = *p.lui_;
    u.t = u.le + 1;
    std::cerr << "t=" << u.t << std::endl;
  }
}

void Suspender::Sem::LoopEnd() {
  U& u = *p.lui_;
  if (u.le < u.lb && // le not initialized for current loop
      u.lb <= u.t) { // target is within the loop
    u.le = u.c;
    std::cerr << "le=" << u.lb << std::endl;
  }
  if ((*this)()) {
    u.t = u.lb + 1;
  }
}


Suspender::Suspender() 
  : lu_(1, U(-1,-1)), lui_(lu_.begin()), nest_(false)
{}

Suspender::Sem Suspender::GetSem(std::string name) {
  return Sem(*this, name);
}

std::string Suspender::GetCurName() const {
  return curname_;
}

std::string Suspender::Print() const {
  std::stringstream b;
  for (auto e : lu_) {
    b << "(" << e.c << " " << e.t << ") ";
  }
  return b.str();
}

bool Suspender::Pending() const {
  return lu_.size() != 1;
}

