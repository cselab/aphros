#pragma once

#include <iostream>
#include <string>
#include <sstream>

#include "Vars.h"

class Interp {
 public:
  Interp() = default;
  void Cmd(std::string s) {
    std::stringstream b(s);
    std::string cmd; 
    b >> cmd;
    if (cmd == "set") {
      CmdSet(s);
    } else {
      std::cerr << "Unknown command: '" << cmd << "'" <<  std::endl;
      assert(false);
    }
  }
  void CmdSet(std::string s) {
    std::string cmd, type, key, val;
    std::stringstream b(s);
    b >> cmd >> type >> key >> val;
    v_.Parse(val, type, key);
  }
  bool Next(std::istream& in) {
    std::string s;
    std::getline(in, s);
    if (!s.empty()) {
      Cmd(s);
      return true;
    }
    return false;
  }
  void All(std::istream& in) {
    while (Next(in)) {}
  }
  template <class T>
  void Print(Vars::Map<T>& m) {
    for (auto& a : m) {
      std::cout 
          << m.GetTypeName() << " " 
          <<  a.first << " = " 
          << m.Print(a.first) << std::endl;
    }
  }
  void PrintAll() {
    Print(v_.String);
    Print(v_.Int);
    Print(v_.Double);
    Print(v_.Vect);
  }

 private:
  Vars v_;
};
