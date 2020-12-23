// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "idx.h"

template <size_t dim_>
class GDir {
 public:
  using MIdx = GMIdx<dim_>;
  static constexpr size_t dim = dim_;

  constexpr GDir() {}
  constexpr explicit GDir(size_t dir) : dir_(dir) {}
  constexpr char GetLetter() const {
    return "xyzw"[dir_];
  }
  explicit operator size_t() const {
    return dir_;
  }
  explicit operator MIdx() const {
    return MIdx::GetUnit(dir_);
  }
  bool operator==(const GDir& other) const {
    return dir_ == other.dir_;
  }
  bool operator!=(const GDir& other) const {
    return !((*this) == other);
  }
  bool operator<(const GDir& other) const {
    return dir_ < other.dir_;
  }
 private:
  size_t dir_;
};
