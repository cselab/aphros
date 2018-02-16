#include <cassert>

#include "suspender.h"

Suspender::Sem::Sem(Suspender& p, std::string name) 
  : p(p), name_(name)
{
  auto& l = p.lu_;
  auto& i = p.lui_;
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

bool Suspender::Sem::operator()() {
  auto& i = p.lui_;
  return i->c++ == i->t;
}

Suspender::Suspender() 
  : lu_(1, U(-1,-1)), lui_(lu_.begin()) 
{}

Suspender::Sem Suspender::GetSem(std::string name) {
  return Sem(*this, name);
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

