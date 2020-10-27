// Created by Petr Karnakov on 07.08.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <iosfwd>
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
  // If `path` is a relative path, looks for the file relative to:
  // - current directory
  // - directory `dir` (if `dir != ""`)
  // Reports file path and line number in case of error.
  void ParseFile(std::string path, std::string dir = "");
  // Prints commands "set" for variables in `map`.
  template <class T>
  static void PrintMap(const Vars::Map<T>& map, std::ostream& out);
  // Prints commands "set" for variables in `var`.
  static void PrintVars(const Vars& var, std::ostream& out);
  // Prints commands "set" for variables in all maps.
  void PrintVars(std::ostream& out) const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
