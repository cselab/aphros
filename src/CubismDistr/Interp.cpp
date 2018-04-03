#include <sstream>
#include <cassert>
#include <fstream>
#include <stdexcept>

#include "Interp.h"
#include "Vars.h"

Interp::Interp(Vars& v) : v_(v) {}

std::string Interp::RemoveComment(std::string s) {
  return s.substr(0, s.find('#', 0));
}

void Interp::Cmd(std::string s) {
  s = RemoveComment(s);
  std::stringstream b(s);
  std::string cmd; 
  b >> cmd;
  if (cmd == "set") {
    CmdSet(s);
  } else if (cmd == "setnext") {
    CmdSetNext(s);
  } else if (cmd == "del") {
    CmdDel(s);
  } else if (cmd == "include") {
    CmdInclude(s);
  } else if (cmd == "") {
    // nop
  } else {
    throw std::runtime_error("Cmd(): unknown command '" + cmd + "'");
  }
}

void Interp::CmdSet(std::string s) {
  std::string cmd, type, key, val;
  std::stringstream b(s);
  b >> std::skipws;
  b >> cmd >> type >> key;
  char c;
  // Read first non-ws character 
  b >> c;
  // Read remaining line
  std::getline(b, val);
  val = c + val;

  v_.SetStr(type, key, val);
}

void Interp::CmdSetNext(std::string s) {
  std::string cmd, type, key, val;
  std::stringstream b(s);
  b >> std::skipws;
  b >> cmd >> type >> key;
  char c;
  // Read first non-ws character 
  b >> c;
  // Read remaining string
  std::getline(b, val);
  val = c + val;

  if (!v_.Int.Exists(key)) {
    v_.Int.Set(key, 0);
  }
  int& id = v_.Int[key];
  v_.SetStr(type, key + std::to_string(id), val);
  ++id;
}


void Interp::CmdDel(std::string s) {
  std::string cmd, key;
  std::stringstream b(s);
  b >> cmd >> key;
  if (!v_.Del(key)) {
    throw std::runtime_error("CmdDel(): unknown variable '" + key + "'");
  }
}

void Interp::CmdInclude(std::string s) {
  std::string cmd, fn;
  std::stringstream b(s);
  b >> cmd >> fn;

  // Read all content
  std::ifstream f(fn);
  if (!f.good()) {
    throw std::runtime_error("CmdInclude: Can't open '" + fn + "'");
  }
  std::stringstream r;
  r << f.rdbuf();
  f.close();

  RunAll(r);
}

void Interp::Run(std::string s) {
  Cmd(s);
}

void Interp::RunNext(std::istream& in) {
  std::string s;
  std::getline(in, s);
  Cmd(s);
}

void Interp::RunAll(std::istream& in) {
  while (in) {
    RunNext(in);
  }
}

template <class T>
void Interp::Print(Vars::Map<T>& m, std::ostream& out) {
  for (auto& a : m) {
    out 
        << "set "
        << m.GetTypeName() << " " 
        << a.first << " "
        << m.GetStr(a.first) << std::endl;
  }
}

void Interp::PrintAll(std::ostream& out) {
  Print(v_.String, out);
  Print(v_.Int, out);
  Print(v_.Double, out);
  Print(v_.Vect, out);
}

void Interp::PrintAll() {
  PrintAll(std::cout);
}

