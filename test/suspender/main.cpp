#include <cassert>
#include <iostream>

#include "suspender.h"

using S = Suspender;
std::string b; 
int n;    

void E(S& s) {
  auto e = s.GetSem();
  for (int i = 0; i < n; ++i) {
    if (e()) {
      b += "E1";
    }
    if (e()) {
      b += "E2";
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
    b += "|";
  } while (s.Pending());

  std::cerr
      << "'" << b << "' == '" << p << "'" << std::endl;
  assert(b == p);

}


int main() {
  Test();
}
