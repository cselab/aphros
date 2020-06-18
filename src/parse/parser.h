// Created by Petr Karnakov on 07.08.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include "vars.h"

class Parser {
 public:
  Parser(Vars& v);
  Parser(const Parser&) = delete;
  Parser(Parser&&) = delete;
  ~Parser();
  Parser& operator=(const Parser&) = delete;
  Parser& operator=(Parser&&) = delete;
  // Executes single command
  void Run(std::string);
  // Executes single line from stream
  void RunNext(std::istream&);
  // Executes all lines from stream
  void RunAll(std::istream&);
  // Executes all lines from file.
  // Includes line number to error messages.
  void RunAll(std::string path);
  // Prints content of Map
  template <class T>
  static void Print(const Vars::Map<T>& m, std::ostream& out);
  // Prints content of v_
  void PrintAll(std::ostream& out) const;
  // Prints content of v_ to std::cout
  void PrintAll() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
