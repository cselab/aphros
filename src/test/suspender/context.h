namespace context {

struct Context {
  Context() {
    std::cerr << "Context created" << std::endl;
  }
  ~Context() {
    std::cerr << "Context destroyed" << std::endl;
  }
};

void B(Suspender& s) {
  Context ctx;
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


void Test() {
  Suspender s;

  do {
    A(s);
  } while (s.Pending());
}

} // namespace context
