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

#include "solver/vof.h"
#include "func/init_u.h"
#include "overlap/overlap.h"
#include "young/young.h"

using Scal = double;
using Vect = GVect<Scal, 3>;
using R = Reconst<Scal>;

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

void TestYoung() {
  std::cerr << "young solution" << std::endl;

  YoungParam q;
  young_set(&q);
  young_ini(q);

  Vect x(0);
  Scal pvy, pvz, p, T;
  young_fields(x[1], x[2], &pvy, &pvz, &p, &T);
  std::cout << std::vector<Scal>{x[1], x[2], pvy, pvz, p, T} << std::endl;
}

int main() {
  TestYoung();
}
