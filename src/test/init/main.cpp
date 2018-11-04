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
#include "overlap/over.h"

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

}

int main() {
  TestYoung();
}
