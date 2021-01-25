// Created by Petr Karnakov on 17.10.2019
// Copyright 2019 ETH Zurich

namespace loop {

using S = Suspender;
std::string b;

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
    std::cerr << s.Print() << " " << s.GetNameSequence() << std::endl;
  } while (s.Pending());

  std::cerr << "'" << b << "' == '" << p << "'" << std::endl;
  assert(b == p);
}

} // namespace loop
