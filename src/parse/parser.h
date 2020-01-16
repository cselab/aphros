// Created by Petr Karnakov on 07.08.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <iostream>
#include <string>

#include "vars.h"

class Parser {
 public:
  Parser(Vars& v);
  // Executes single command
  void Run(std::string);
  // Executes single line from stream
  void RunNext(std::istream&);
  // Executes all lines from stream
  void RunAll(std::istream&);
  // Prints content of Map
  template <class T>
  static void Print(const Vars::Map<T>& m, std::ostream& out);
  // Prints content of v_
  void PrintAll(std::ostream& out) const;
  // Prints content of v_ to std::cout
  void PrintAll() const;

 private:
  static std::string RemoveComment(std::string s);
  // Executes single command
  void Cmd(std::string);
  // set <type> <key> <value>
  void CmdSet(std::string);
  // setnext <type> <key> <value>
  // Sets value with name key+id and increments id stored in v_.Int[key]
  void CmdSetNext(std::string);
  // del <name>
  void CmdDel(std::string);
  // include <filename>
  void CmdInclude(std::string);

  Vars& v_;
};
