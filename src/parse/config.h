// Created by Petr Karnakov on 01.05.2020
// Copyright 2020 ETH Zurich

#pragma once

#include "geom/vect.h"
#include "vars.h"

#define VAR_GENERIC(TYPE, PTRS, NAME) \
  TYPE NAME;                          \
  void* NAME##_reg_ = (PTRS[#NAME] = &NAME)

#define VAR_DOUBLE(NAME) VAR_GENERIC(double, ptrs_double_, NAME)
#define VAR_INT(NAME) VAR_GENERIC(int, ptrs_int_, NAME)
#define VAR_STRING(NAME) VAR_GENERIC(std::string, ptrs_string_, NAME)
#define VAR_VECT(NAME) VAR_GENERIC(std::vector<double>, ptrs_vect_, NAME)
#define VAR_BOOL(NAME) VAR_GENERIC(bool, ptrs_bool_, NAME)
#define VAR_VECT3(NAME) VAR_GENERIC(Vect3, ptrs_vect3_, NAME)

class ConfigBase {
 public:
  ConfigBase() = default;
  ConfigBase(const ConfigBase&) = delete;
  ConfigBase(ConfigBase&&) = delete;
  ConfigBase& operator=(const ConfigBase&) = delete;
  ConfigBase& operator=(ConfigBase&&) = delete;
  void Read(const Vars& var, const std::string prefix = "") {
    for (auto p : ptrs_double_) {
      *p.second = var.Double[prefix + p.first];
    }
    for (auto p : ptrs_int_) {
      *p.second = var.Int[prefix + p.first];
    }
    for (auto p : ptrs_string_) {
      *p.second = var.String[prefix + p.first];
    }
    for (auto p : ptrs_vect_) {
      *p.second = var.Vect[prefix + p.first];
    }
    for (auto p : ptrs_bool_) {
      *p.second = var.Int[prefix + p.first];
    }
    for (auto p : ptrs_vect3_) {
      *p.second = Vect3(var.Vect[prefix + p.first]);
    }
  }

 protected:
  using Vect3 = generic::Vect<double, 3>;
  std::map<std::string, double*> ptrs_double_;
  std::map<std::string, int*> ptrs_int_;
  std::map<std::string, std::string*> ptrs_string_;
  std::map<std::string, std::vector<double>*> ptrs_vect_;
  std::map<std::string, bool*> ptrs_bool_;
  std::map<std::string, Vect3*> ptrs_vect3_;
};
