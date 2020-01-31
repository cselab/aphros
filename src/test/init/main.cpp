// Created by Petr Karnakov on 04.11.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "func/init_u.h"
#include "overlap/overlap.h"
#include "solver/vofm.h"
#include "young/young.h"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;
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
