#undef NDEBUG
#include <iostream>
#include <string>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <functional>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <random>

#include "util/height.h"
#include "debug/isnan.h"

using Scal = double;
using U = UHeight<Scal>;

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

void TestGood() {
  std::cout << "UHeight::Good(v,n)" << std::endl;
  auto t = [](const std::vector<Scal>& v, Scal n, bool ref) {
    auto b = (!IsNan(U::Good(v, n)) ? "true" : "false");
    auto br = (ref ? "true" : "false");
    auto eq = (b == br ? "==" : "!=");
    std::cout << v << " n=" << n;
    std::cout << " , " << b << eq << br << std::endl;
    assert(b == br);
  };
  t({0, 0, 0, 0.5, 1, 1, 1}, -1, true);
  t({1, 1, 1, 0.5, 0, 0, 0}, 1, true);
  t({1, 0.5, 0.5, 0.5, 0.5, 0.5, 0}, 1, true);
  t({1, 0.7, 0.6, 0.5, 0.4, 0.3, 0}, 1, true);
  t({0, 0.7, 0.6, 0.5, 0.4, 0.3, 1}, 1, false);
  t({1, 1, 1, 0.5, 0, 0, 0}, 1, true);
  t({0, 0, 0, 0.5, 1, 1, 1}, 1, false);
  t({1, 0.5, 0, 0, 0, 0, 0}, 1, true);
  t({1, 0.5, 1, 1, 0, 0, 0}, 1, false);
  t({0, 0.5, 0, 0, 1, 1, 1}, 1, false);
  t({0, 0.5, 1, 0.5, 0, 0, 0}, 1, true);
  t({0, 0, 0, 0.5, 1, 0.5, 0}, 1, false);
  t({0, 0, 0, 0.5, 1, 0.5, 0}, -1, true);
  t({1, 1, 1, 0.5, 0.8, 0, 0}, 1, true);
  t({1,1,0.42,0,0,0,0}, 1, true);
}

void TestGood2() {
  std::cout << "UHeight::Good(v)" << std::endl;
  auto t = [](const std::vector<Scal>& v, Scal ref) {
    Scal r = U::Good(v);
    bool cl = (std::abs(r - ref) < 1e-6);
    if (IsNan(r)) {
      cl = IsNan(ref);
    }
    auto eq = (cl ? "==" : "!=");
    std::cout << v;
    std::cout << " , " << r << eq << ref << std::endl;
    assert(cl);
  };
  Scal nan = GetNan<Scal>();
  t({0, 0, 0, 0.5, 1, 1, 1}, 0);
  t({1, 1, 1, 0.5, 0, 0, 0}, 0);
  t({1, 0.5, 0.5, 0.5, 0.5, 0.5, 0}, 0);
  t({1, 0.7, 0.6, 0.5, 0.4, 0.3, 0}, 0);
  t({0, 0.7, 0.6, 0.5, 0.4, 0.3, 1}, 0);
  t({1, 1, 1, 0.1, 0, 0, 1}, -0.4);
  t({0, 0, 0, 0.1, 1, 1, 1}, 0.4);
  t({1, 0.5, 0.5, 1, 0.5, 0.5, 0}, nan);
  t({0, 0.5, 0, 0, 1, 1, 1}, nan);
  t({0, 0.5, 0, 0.1, 1, 1, 1}, 0.4);
  t({0, 0, 0, 0.5, 1, 0.5, 0}, 0);
  t({1, 1, 1, 0.5, 0.8, 0, 0}, 0.8);
  t({1,1,0.42,0.1,0,0,0}, -0.98);
   
}


int main() {
  TestGood();
  TestGood2();
}
