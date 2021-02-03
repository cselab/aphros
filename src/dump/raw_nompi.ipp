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
void Raw<M>::WriteXmf(
    std::string, std::string, Type, Vect, Vect, MIdx, std::string) {}

template <class M>
void Raw<M>::WriteXmf(std::string, std::string, Type, std::string, const M&) {}

} // namespace dump
