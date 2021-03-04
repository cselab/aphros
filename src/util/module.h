// Created by Petr Karnakov on 20.04.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include "util/logger.h"

#define DECLARE_FORCE_LINK_TARGET(x) int kForceLink_##x = 0
#define FORCE_LINK(x)          \
  do {                         \
    extern int kForceLink_##x; \
    kForceLink_##x = 1;        \
  } while (0)

template <class Mod>
bool RegisterModule() {
  Mod* ptr = new Mod();
  const auto name = ptr->GetName();
  fassert(
      !Mod::GetInstance(name),
      "RegisterModule: module '" + name + "' of type '" + typeid(Mod).name() +
          "' and base type '" + typeid(typename Mod::Base).name() +
          "' already registered");
  Mod::AddInstance(name, ptr);
  return true;
}

// XXX Note about namespaces
// Definitions of module implementations must be inside a separate namespace
// for each module type.
// Otherwise, calls such as `RegisterModule<Uniform<M>>()`
// instantiate functions that do not contain the module type in their
// signature.
// Therefore, two classes implementing modules of different types
// (e.g. ModuleInitVelocity and ModuleInitContang)
// but having the same name (e.g. `Uniform`)
// would generate symbols of the same name and cause a linker symbol clash.

template <class Base_>
class Module {
 public:
  using Base = Base_;
  Module(std::string name, std::string desc = "") : name_(name), desc_(desc) {}
  virtual ~Module() = default;
  std::string GetName() const {
    return name_;
  }
  std::string GetDesc() const {
    return desc_;
  }
  static std::map<std::string, Base*> GetInstances() {
    std::map<std::string, Base*> map;
    for (auto& p : GetTable()) {
      map[p.first] = p.second.get();
    }
    return map;
  }
  static void AddInstance(std::string name, Base* ptr) {
    GetTable()[name].reset(ptr);
  }
  static Base* GetInstance(std::string name) {
    auto& table = GetTable();
    auto it = table.find(name);
    if (it != table.end()) {
      return it->second.get();
    }
    return nullptr;
  }

 private:
  const std::string name_;
  const std::string desc_;
  // Access through GetTable() ensures that the table is constructed on first
  // use (would be undefined for member variables).
  static auto& GetTable() {
    static std::map<std::string, std::unique_ptr<Base>> table;
    return table;
  }
};
