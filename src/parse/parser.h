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
  // Returns command `set TYPE KEY VALUE`
  // where `KEY` is `key`, and `TYPE` and `VALUE` correspond
  // to a variable `key` found in `var`.
  // If `new_key` is not empty, replaces `KEY` with `new_key.`
  // If `key` is not found, throws exception.
  static std::string GenerateSetCommand(
      const Vars& var, const std::string& key, const std::string new_key = "");
  // Writes a new line with command returned by `GenerateSetCommand()` to stream.
  static void WriteSetCommand(
      std::ostream& out, const Vars& var, const std::string& key,
      const std::string new_key = "");

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
