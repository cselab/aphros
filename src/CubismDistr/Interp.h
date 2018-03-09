#pragma once

#include <iostream>
#include <string>

#include "Vars.h"

class Interp {
 public:
  Interp(Vars& v);
  // Executes single line 
  void RunNext(std::istream& in); 
  // Executes all lines 
  void RunAll(std::istream& in); 
  // Prints content of Map
  template <class T>
  static void Print(Vars::Map<T>& m, std::ostream& out); 
  // Prints content of v_
  void PrintAll(std::ostream& out); 
  // Prints content of v_ to std::cout
  void PrintAll(); 

 private:
  static std::string RemoveComment(std::string s);
  // Executes single command
  void Cmd(std::string s);
  // set <type> <key> <value>
  void CmdSet(std::string s);
  // del <name>
  void CmdDel(std::string s); 

  Vars& v_;
};
