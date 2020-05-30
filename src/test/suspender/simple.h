// Created by Petr Karnakov on 17.10.2019
// Copyright 2019 ETH Zurich

namespace simple {

void B(Suspender& s) {
  Suspender::Sem sem = s.GetSem();
  if (sem()) {
    std::cerr << "B1" << std::endl;
  }
  if (sem()) {
    std::cerr << "B2" << std::endl;
  }
}

void A(Suspender& s) {
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
