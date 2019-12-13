#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>
#include <typeinfo>

#include "solver/cond.h"
#include "util/hydro.ipp"

const int dim = 3;
using Scal = double;
using Vect = GVect<Scal, dim>;

using namespace solver;

std::string P(const void*) {
  return std::string();
}

std::string P(const CondFace* b) {
  std::stringstream ss;
  ss << " " << b->GetNci();
  return ss.str();
}

template <class T>
std::string P(const CondFaceVal<T>* d) {
  std::stringstream ss;
  ss << " " << d->GetNci() << " " << d->GetValue();
  return ss.str();
}

template <class T>
std::string P(const CondFaceGrad<T>* d) {
  std::stringstream ss;
  ss << " " << d->GetNci() << " " << d->GetGrad();
  return ss.str();
}

template <class D>
void Try(const UniquePtr<CondFace>& b, std::string& s) {
  if (auto d = b.Get<D>()) {
    s = typeid(d).name() + P(d);
  }
}

template <class V>
class CondFaceValCustom : public CondFaceVal<V> {
 public:
  CondFaceValCustom(const V& v, size_t nci) : CondFaceVal<V>(nci), v_(v) {}
  V GetValue() const override {
    return v_;
  }

 private:
  V v_;
};

template <class V>
class CondFaceGradCustom : public CondFaceGrad<V> {
 public:
  CondFaceGradCustom(const V& v, size_t nci) : CondFaceGrad<V>(nci), v_(v) {}
  V GetGrad() const override {
    return v_;
  }

 private:
  V v_;
};

class CondFaceReflectCustom : public CondFaceReflect {
 public:
  CondFaceReflectCustom(size_t nci) : CondFaceReflect(nci) {}
};

class CondFaceExtrapCustom : public CondFaceExtrap {
 public:
  CondFaceExtrapCustom(size_t nci) : CondFaceExtrap(nci) {}
};

void Print(const UniquePtr<CondFace>& b) {
  std::string s{"none"};
  Try<CondFace>(b, s);
  Try<CondFaceVal<Scal>>(b, s);
  Try<CondFaceVal<Vect>>(b, s);
  Try<CondFaceValFixed<Scal>>(b, s);
  Try<CondFaceValFixed<Vect>>(b, s);
  Try<CondFaceValCustom<Scal>>(b, s);
  Try<CondFaceValCustom<Vect>>(b, s);

  Try<CondFaceGrad<Scal>>(b, s);
  Try<CondFaceGrad<Vect>>(b, s);
  Try<CondFaceGradFixed<Scal>>(b, s);
  Try<CondFaceGradFixed<Vect>>(b, s);
  Try<CondFaceGradCustom<Scal>>(b, s);
  Try<CondFaceGradCustom<Vect>>(b, s);

  Try<CondFaceReflect>(b, s);
  Try<CondFaceExtrap>(b, s);

  std::cout << s << std::endl;
}

#define EPrint(x)                           \
  {                                         \
    std::cout << "-> " << #x << std::endl;  \
    try {                                   \
      Print(x);                             \
    } catch (const std::runtime_error& e) { \
      std::cout << e.what() << std::endl;   \
    }                                       \
  }

#define ECHO(x)                 \
  std::cout << #x << std::endl; \
  x;

void Test() {
  std::cout << "\nTest Eval" << std::endl;
  UniquePtr<CondFace> b;
  ECHO()
  ECHO(b.Set<CondFaceValCustom<Scal>>(3.14, 1);)
  EPrint(b);
  EPrint(Eval<Scal>(b));
  EPrint(Eval<Vect>(b));
  EPrint(EvalComp<Vect>(b, 0));

  ECHO()
  ECHO(b.Set<CondFaceGradCustom<Scal>>(3.14, 1);)
  EPrint(b);
  EPrint(Eval<Scal>(b));
  EPrint(Eval<Vect>(b));
  EPrint(EvalComp<Vect>(b, 0));

  ECHO()
  ECHO(b.Set<CondFaceValCustom<Vect>>(Vect(1.1, 2.2, 3.3), 1);)
  EPrint(b);
  EPrint(Eval<Scal>(b));
  EPrint(Eval<Vect>(b));
  EPrint(EvalComp<Vect>(b, 0));
  EPrint(EvalComp<Vect>(b, 1));
  EPrint(EvalComp<Vect>(b, 2));

  ECHO()
  ECHO(b.Set<CondFaceGradCustom<Vect>>(Vect(1.1, 2.2, 3.3), 1);)
  EPrint(b);
  EPrint(Eval<Scal>(b));
  EPrint(Eval<Vect>(b));
  EPrint(EvalComp<Vect>(b, 0));
  EPrint(EvalComp<Vect>(b, 1));
  EPrint(EvalComp<Vect>(b, 2));

  ECHO()
  ECHO(b.Set<CondFaceReflectCustom>(1);)
  EPrint(b);
  EPrint(Eval<Scal>(b));
  EPrint(Eval<Vect>(b));

  ECHO()
  ECHO(b.Set<CondFaceExtrapCustom>(1);)
  EPrint(b);
  EPrint(Eval<Scal>(b));
  EPrint(Eval<Vect>(b));
}

void TestParse() {
  std::cout << "\nTest Parse" << std::endl;
  CondFaceAdvection<Scal> cfa;
  cfa.nci = 1;
  ParseAdvectionFaceCond("clear0 1", cfa);
  ParseAdvectionFaceCond("clear1 2", cfa);
  ParseAdvectionFaceCond("halo fill", cfa);
  ParseAdvectionFaceCond("fill_vf 0.4", cfa);
  ParseAdvectionFaceCond("fill_cl 0.6", cfa);
  std::cout << cfa << std::endl;
  ParseAdvectionFaceCond("halo reflect", cfa);
  std::cout << cfa << std::endl;
  for (auto s : Split("clear0 2, clear1 3, fill_vf 0.3, fill_cl 0.7", ',')) {
    ParseAdvectionFaceCond(s, cfa);
  }
  std::cout << cfa << std::endl;
}

int main() {
  Test();
  TestParse();
}
