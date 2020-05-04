// Created by Petr Karnakov on 20.04.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

template <class Mod>
bool RegisterModule() {
  Mod* ptr = new Mod();
  const auto name = ptr->GetName();
  if (Mod::GetInstance(name)) {
    throw std::runtime_error(
        "RegisterModule: module '" + name + "' already registered");
  }
  Mod::AddInstance(name, ptr);
  return true;
}

template <class Base>
class Module {
 public:
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
  // Access through GetTable() ensures that table is constructed on first use
  // (in contrast to a private member variable).
  static auto& GetTable() {
    static std::map<std::string, std::unique_ptr<Base>> table;
    return table;
  }
};
