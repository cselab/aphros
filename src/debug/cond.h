// Created by Petr Karnakov on 21.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <sstream>
#include <string>
#include <typeinfo>

#include "solver/cond.h"

namespace debugimpl {

const int dim = 3;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;

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
  ss << " " << d->GetNci() << " " << d->second();
  return ss.str();
}

template <class T>
std::string P(const CondFaceGrad<T>* d) {
  std::stringstream ss;
  ss << " " << d->GetNci() << " " << d->GetGrad();
  return ss.str();
}

template <class D>
void Try(const CondFace* b, std::string& s) {
  if (auto d = dynamic_cast<const D*>(b)) {
    s = typeid(d).name() + P(d);
  }
}

std::string GetName(const CondFace* b) {
  std::string s{"none"};
  Try<CondFace>(b, s);
  Try<CondFaceVal<Scal>>(b, s);
  Try<CondFaceVal<Vect>>(b, s);
  Try<CondFaceValFixed<Scal>>(b, s);
  Try<CondFaceValFixed<Vect>>(b, s);

  Try<CondFaceGrad<Scal>>(b, s);
  Try<CondFaceGrad<Vect>>(b, s);
  Try<CondFaceGradFixed<Scal>>(b, s);
  Try<CondFaceGradFixed<Vect>>(b, s);

  Try<CondFaceReflect>(b, s);
  Try<CondFaceExtrap>(b, s);
  return s;
}

} // namespace debugimpl

std::string GetCondName(const CondFace* b) {
  return debugimpl::GetName(b);
}
