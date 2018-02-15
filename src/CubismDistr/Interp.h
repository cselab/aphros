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
    } else if (cmd == "del") {
      CmdDel(s);
    } else if (cmd == "") {
      // nop
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
  template <class T>
  bool Del(Vars::Map<T>& m, std::string k) {
    if (m.Exists(k)) {
      m.Del(k);
      return true;
    }
    return false;
  }
  bool Del(std::string k) {
    bool r = false;
    if (!r) r = Del(v_.String, k);
    if (!r) r = Del(v_.Int, k);
    if (!r) r = Del(v_.Double, k);
    if (!r) r = Del(v_.Vect, k);
    return r;
  }
  void CmdDel(std::string s) {
    std::string cmd, key;
    std::stringstream b(s);
    b >> cmd >> key;
    if (!Del(key)) {
      std::cerr << "del: unknown variable '" << key << "'" << std::endl;
      assert(false);
    }
  }
  bool Next(std::istream& in) {
    std::string s;
    std::getline(in, s);
    if (in) {
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
