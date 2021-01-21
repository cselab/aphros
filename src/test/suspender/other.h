// Created by Petr Karnakov on 17.10.2019
// Copyright 2019 ETH Zurich

namespace other {

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

void C(S&) { // no sem
  b += "C1";
  b += "C2";
}

void BB(S& s) {
  auto e = s.GetSem("BB"); // dummy sem
  e();
}

void B(S& s) {
  // s.GetSem();  // dummy sem
  auto e = s.GetSem("B");
  // s.GetSem();  // dummy sem
  if (e("1")) {
    b += "B1";
  }
  if (e.Nested("C")) {
    C(s);
  }
  // s.GetSem();  // dummy sem in the middle
  if (e("2")) {
    b += "B2";
  }
}

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
    std::cerr << s.Print() << " " << s.GetNameSequence() << std::endl;
  } while (s.Pending());

  std::cerr << "'" << b << "' == '" << p << "'" << std::endl;
  assert(b == p);
}

} // namespace other
