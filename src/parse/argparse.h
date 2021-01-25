// Created by Petr Karnakov on 26.09.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <cctype>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "parse/vars.h"
#include "util/format.h"
#include "util/logger.h"

class ArgumentParser {
 public:
  template <class T>
  class Proxy {
   public:
    Proxy() = delete;
    Proxy(ArgumentParser& parser, std::string key);
    Proxy& Help(std::string help);
    Proxy& Options(const std::vector<T>&);

   private:
    ArgumentParser& parser_;
    const std::string key_;
  };

  ArgumentParser(std::string desc = "", bool isroot = true);
  ~ArgumentParser();
  // Adds a swtich.
  // If present in arguments, known_args_.Int[key] is set to 1, otherwise 0.
  Proxy<int> AddSwitch(const std::vector<std::string>& names);
  Proxy<int> AddSwitch(std::initializer_list<std::string> names) {
    return AddSwitch(std::vector<std::string>{names});
  }
  Proxy<int> AddSwitch(std::string name) {
    return AddSwitch({name});
  }
  template <class T>
  Proxy<T> AddVariable(
      const std::vector<std::string>& names, T defaultval, bool hasdefault);
  template <class T>
  auto AddVariable(std::initializer_list<std::string> names, T defaultval) {
    return AddVariable<T>(std::vector<std::string>{names}, defaultval, true);
  }
  template <class T>
  auto AddVariable(std::string name, T defaultval) {
    return AddVariable<T>({name}, defaultval, true);
  }
  template <class T>
  auto AddVariable(std::initializer_list<std::string> names) {
    return AddVariable<T>(names, {}, false);
  }
  template <class T>
  auto AddVariable(std::string name) {
    return AddVariable<T>({name}, {}, false);
  }
  const Vars& GetKnownArgs() const;
  Vars ParseArgs(std::vector<std::string> argv, std::string program = "") const;
  Vars ParseArgs(
      int argc, const char** argv, std::string ignore_after = "") const;
  void PrintHelp(std::ostream& out, bool full, std::string program) const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
