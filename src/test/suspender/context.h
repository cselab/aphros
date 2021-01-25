// Created by Petr Karnakov on 17.10.2019
// Copyright 2019 ETH Zurich

#include <vector>

namespace context {

void B(Suspender& susp) {
  struct Context {
    Context() {
      std::cerr << "ContextB created" << std::endl;
    }
    ~Context() {
      std::cerr << "ContextB destroyed" << std::endl;
    }
    int a = 10;
  };
  Suspender::Sem sem = susp.GetSem();
  Context* ctx(sem);
  if (sem()) {
    ++ctx->a;
    std::cerr << "B1" << std::endl;
  }
  if (sem()) {
    ++ctx->a;
    std::cerr << "B2 " << ctx->a << std::endl;
  }
}

void A(Suspender& susp) {
  struct Context {
    Context() {
      std::cerr << "ContextA created" << std::endl;
    }
    ~Context() {
      std::cerr << "ContextA destroyed" << std::endl;
    }
    std::vector<double> v;
  };
  Suspender::Sem sem = susp.GetSem();
  Context* ctx(sem);
  if (sem()) {
    ctx->v = {1., 2., 3., 4.};
    std::cerr << "A1" << std::endl;
  }
  if (sem.Nested()) {
    B(susp);
  }
  if (sem.Nested()) {
    B(susp);
  }
  if (sem()) {
    std::cerr << "A2 " << ctx->v.size() << std::endl;
  }
}

void Test() {
  Suspender susp;

  do {
    A(susp);
  } while (susp.Pending());
}

} // namespace context
