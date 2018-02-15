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

#include <iostream>

namespace test_suspender {

using S = Suspender;
std::string b; // init in Test()
int n;         // init in Test()

void E(S& s) {
  auto e = s.GetSem();
  for (int i = 0; i < n; ++i) {
    if (e()) {
      b += "E1";
    }
    if (e()) {
      b += "E2";
      --n;
    }
  }
}

void D(S& s) {
  auto e = s.GetSem(); // sem but no stages
  b += "D1";
  b += "D2";
}

void C(S& s) { // no sem
  b += "C1";
  b += "C2";
}

void B(S& s) {
  auto e = s.GetSem();
  s.GetSem();  // dummy sem
  if (e()) {
    b += "B1";
  }
  if (e()) {
    C(s);
  }
  s.GetSem();  // dummy sem in the middle
  if (e()) {
    b += "B2";
  }
}

void A(S& s) {
  auto e = s.GetSem();
  if (e()) {
    b += "A1";
  }
  if (e()) {
    B(s);
  }
  if (e()) {
    b += "A2";
  }
  if (e()) {
    C(s);
  }
  if (e()) {
    D(s);
  }
  if (e()) {
    E(s);
  }
}

void Test() {
  S s;
  b = "";
  n = 3;
  std::string p = "A1|B1|C1C2|B2|A2|C1C2|D1D2|E1|E2|E1|E2|E1|E2|";
  do {
    A(s);
    std::cerr << s.Print() << std::endl;
    b += "|";
  } while (s.Pending());

  std::cerr
      << "'" << b << "' == '" << p << "'" << std::endl;
  assert(b == p);

}

} // namespace test_suspender
