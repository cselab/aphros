// Created by Petr Karnakov on 01.02.2021
// Copyright 2021 ETH Zurich

#include <fstream>
#include <limits>
#include <sstream>

#include "raw.h"
#include "util/distr.h"
#include "util/format.h"
#include "util/logger.h"

namespace dump {

template <class M>
template <class T>
void Raw<M>::Write(const FieldCell<T>&, const Meta&, std::string, M&) {}

template <class M>
template <class T>
void Raw<M>::Read(FieldCell<T>&, const Meta&, std::string, M&) {}

template <class M>
std::string GetXmfTemplate();

template <class M>
auto Raw<M>::GetMeta(MIdx start, MIdx stride, const M&) -> Meta {
  return {};
}
template <class M>
void Raw<M>::WriteXmf(std::ostream&, const Meta&) {}
template <class M>
void Raw<M>::WriteXmf(const std::string& xmfpath, const Meta&) {}

template <class M>
auto Raw<M>::ReadXmf(std::istream&) -> Meta {
  return {};
}
template <class M>
auto Raw<M>::ReadXmf(const std::string& xmfpath) -> Meta {
  return {};
}

} // namespace dump
