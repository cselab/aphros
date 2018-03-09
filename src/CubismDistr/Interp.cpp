#include <sstream>
#include <cassert>

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
  } else if (cmd == "del") {
    CmdDel(s);
  } else if (cmd == "") {
    // nop
  } else {
    std::cerr << "Unknown command: '" << cmd << "'" <<  std::endl;
    assert(false);
  }
}

void Interp::CmdSet(std::string s) {
  std::string cmd, type, key, val;
  std::stringstream b(s);
  b >> std::skipws;
  b >> cmd >> type >> key;
  // Read first non-ws character (append later)
  char c;
  b >> c;
  std::getline(b, val);
  v_.Parse(c + val, type, key);
}

template <class T>
bool Interp::Del(Vars::Map<T>& m, std::string k) {
  if (m.Exists(k)) {
    m.Del(k);
    return true;
  }
  return false;
}

bool Interp::Del(Vars& v, std::string k) {
  bool r = false;
  if (!r) r = Del(v.String, k);
  if (!r) r = Del(v.Int, k);
  if (!r) r = Del(v.Double, k);
  if (!r) r = Del(v.Vect, k);
  return r;
}

void Interp::CmdDel(std::string s) {
  std::string cmd, key;
  std::stringstream b(s);
  b >> cmd >> key;
  if (!Del(v_, key)) {
    std::cerr << "del: unknown variable '" << key << "'" << std::endl;
    assert(false);
  }
}

bool Interp::RunNext(std::istream& in) {
  std::string s;
  std::getline(in, s);
  if (in) {
    Cmd(s);
    return true;
  }
  return false;
}

void Interp::RunAll(std::istream& in) {
  while (RunNext(in)) {}
}

template <class T>
void Interp::Print(Vars::Map<T>& m, std::ostream& out) {
  for (auto& a : m) {
    out 
        << m.GetTypeName() << " " 
        <<  a.first << " = " 
        << m.Print(a.first) << std::endl;
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

