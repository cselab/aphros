#pragma once

#include <iostream>
#include <string>

#include "Vars.h"

class Interp {
 public:
  Interp(Vars& v);
  std::string RemoveComment(std::string s);
  void Cmd(std::string s);
  void CmdSet(std::string s);
  template <class T>
  static bool Del(Vars::Map<T>& m, std::string k); 
  static bool Del(Vars& v, std::string k); 
  void CmdDel(std::string s); 
  bool RunNext(std::istream& in); 
  void RunAll(std::istream& in); 
  template <class T>
  void Print(Vars::Map<T>& m, std::ostream& out); 
  void PrintAll(std::ostream& out); 
  void PrintAll(); 

 private:
  Vars& v_;
};
