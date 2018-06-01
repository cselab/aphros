#undef NDEBUG
#include <cassert>
#include <iostream>

#include "util/suspender.h"

namespace simple {

void B (Suspender& s) {
  Suspender::Sem sem = s.GetSem();
  if (sem()) {
    std::cerr << "B1" << std::endl;
  }
  if (sem()) {
    std::cerr << "B2" << std::endl;
  }
}

void A (Suspender& s) {
  Suspender::Sem sem = s.GetSem();
  if (sem()) {
    std::cerr << "A1" << std::endl;
  }
  if (sem.Nested()) {
    B(s);
  }
  if (sem()) {
    std::cerr << "A2" << std::endl;
  }
}


void Simple() {
  Suspender s;

  do {
    A(s);
  } while (s.Pending());
}

} // namespace simple


using S = Suspender;

std::string b; 

void E(S& s) {
  auto e = s.GetSem("E");
  for (int i = 0; i < 3; ++i) {
    if (e("1")) {
      b += "E1";
    }
    if (e("2")) {
      b += "E2";
    }
  }
}

void D(S& s) {
  auto e = s.GetSem("D"); // sem but no stages
  b += "D1";
  b += "D2";
}

void C(S& s) { // no sem
  b += "C1";
  b += "C2";
}

void BB(S& s) {
  auto e = s.GetSem("BB");  // dummy sem 
  e();
}

void B(S& s) {
  //s.GetSem();  // dummy sem 
  auto e = s.GetSem("B");
  //s.GetSem();  // dummy sem
  if (e("1")) {
    b += "B1";
  }
  if (e.Nested("C")) {
    C(s);
  }
  //s.GetSem();  // dummy sem in the middle
  if (e("2")) {
    b += "B2";
  }
}

namespace loop {

int k;

// action that changes conv condition
void Iter(S& s) {
  auto e = s.GetSem("I");
  if (e("a")) {
    b += "a";
  }
  if (e("+")) {
    b += "+";
    ++k;
  }
  if (e("b")) {
    b += "b";
  }
  std::cerr << "k=" << k << std::endl;
}

// Loop suspender:
// * syntax of loop must be preserved
// * additional calls in either condition check 
// or begin/end of loop body
// * allow multi-stage loop body
// * start with specific form
//   while (diff > tol) {
//     
//   }
//   - both diff and tol are global variables
//   - diff is modified within the loop body
//


void Run(S& s) {
  auto e = s.GetSem("L");
  if (e()) {
    b += "0";
    k = 0;
  }
  if (e()) {
    b += "lb";
  }
  e.LoopBegin();
  if (e.Nested()) {
    Iter(s);
  }
  if (e()) {
    b += "c";
    if (k >= 3) {
      e.LoopBreak();
    }
  }
  e.LoopEnd();

  if (e()) {
    b += "le";
  }
}

void Test() {
  S s;
  b = "";
  std::string p = "0|lb|a|+|b|c|a|+|b|c|a|+|b|c|le|";
  do {
    Run(s);
    b += "|";
    std::cerr << s.Print() << " " << s.GetCurName() << std::endl;
  } while (s.Pending());

  std::cerr
      << "'" << b << "' == '" << p << "'" << std::endl;
  assert(b == p);
}

} // namespace loop


void A(S& s) {
  auto e = s.GetSem("A");
  if (e("1")) {
    b += "A1";
  }
  if (e.Nested("B")) {
    // BB(s);  // TODO: would cause infinite loop (see todo in suspender.h)
    B(s);
  }
  if (e("2")) {
    b += "A2";
  }
  if (e.Nested("C")) {
    C(s);
  }
  if (e.Nested("D")) {
    D(s);
  }
  if (e.Nested("E")) {
    E(s);
  }
}

void Test() {
  S s;
  b = "";
  std::string p = "A1|B1|C1C2|B2|A2|C1C2|D1D2|E1|E2|E1|E2|E1|E2|";
  do {
    A(s);
    b += "|";
    std::cerr << s.Print() << " " << s.GetCurName() << std::endl;
  } while (s.Pending());

  std::cerr
      << "'" << b << "' == '" << p << "'" << std::endl;
  assert(b == p);

}

int main() {
  simple::Simple();

  Test();
  loop::Test();
}
