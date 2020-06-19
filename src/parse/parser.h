// Created by Petr Karnakov on 07.08.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

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
  void ParseStream(std::istream&);
  // Executes all lines from file.
  // Includes line number to error messages.
  void ParseFile(std::string path);
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

// Parses arguments.
//   argc, argv: arguments
//   novalue: names of arguments without value (which is set to "")
// Returns:
//   - map `key:value` for arguments starting from '--' or '-'
//   - list of positional arguments
// Example:
// ./main --a 0 -b -c 2 a b c
// {
//   {
//     {"--a", "0"},
//     {"-b", ""},
//     {"-c", "2"},
//   },
//   {"a", "b", "c"},
// }
std::pair<std::map<std::string, std::string>, std::vector<std::string>>
ParseArgs(int argc, const char** argv, const std::set<std::string>& novalue);
