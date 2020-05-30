// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <map>

#include "field.h"

template <class Idx>
struct GIdxHash {
  size_t operator()(Idx i) const noexcept {
    return size_t(i);
  }
};

template <class Idx>
struct GIdxLess {
  bool operator()(Idx a, Idx b) const noexcept {
    return size_t(a) < size_t(b);
  }
};

// Idx_: instance of GIdx
template <class Value_, class Idx_>
class GMap {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Container = std::map<Idx, Value, GIdxLess<Idx>>;

  size_t size() const {
    return d_.size();
  }
  void clear() {
    d_.clear();
  }
  Value& operator[](Idx i) {
    return d_[i];
  }
  Value& at(Idx i) {
    return d_.at(i);
  }
  const Value& at(Idx i) const {
    return d_.at(i);
  }
  Value* find(Idx i) {
    auto it = d_.find(i);
    return it != d_.end() ? &it->second : nullptr;
  }
  const Value* find(Idx i) const {
    auto it = d_.find(i);
    return it != d_.end() ? &it->second : nullptr;
  }
  void erase(Idx i) {
    d_.erase(i);
  }

  auto begin() {
    return d_.begin();
  }
  auto end() {
    return d_.end();
  }
  auto begin() const {
    return d_.begin();
  }
  auto end() const {
    return d_.end();
  }
  auto cbegin() const {
    return d_.cbegin();
  }
  auto cend() const {
    return d_.cend();
  }

 private:
  Container d_;
};

template <class T>
using MapCell = GMap<T, IdxCell>;

template <class T>
using MapFace = GMap<T, IdxFace>;

template <class T>
using MapNode = GMap<T, IdxNode>;
