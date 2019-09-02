#pragma once

#include <memory>

template <class M_>
class UVof {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  UVof(M&);
  ~UVof();
  void Dump();

 public:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
