// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#include <fstream>
#include "hdf.h"
#include "util/logger.h"

namespace dump {

template <class M>
template <class Field>
void Hdf<M>::Write(const Field&, std::string, M&, std::string) {}

template <class M>
template <class Field>
void Hdf<M>::Read(Field&, std::string, M&, std::string) {}

template <class M>
std::vector<size_t> Hdf<M>::GetShape(std::string, std::string) {
  return {0, 0, 0};
}

template <class M>
void Hdf<M>::WriteXmf(
    std::string, std::string, typename M::Vect, typename M::Vect,
    typename M::MIdx, std::string, std::string) {}

template <class M>
void Hdf<M>::WriteXmf(
    std::string, std::string, std::string, const M&, std::string) {}

template <class M>
void Hdf<M>::WriteBlocks(
    const std::string&, const std::vector<MIdx>&, const std::vector<MIdx>&,
    const std::vector<std::vector<Scal>>&, MIdx, Type, std::string,
    const MpiWrapper&, bool) {}

} // namespace dump
