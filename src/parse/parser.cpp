// Created by Petr Karnakov on 07.08.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "parser.h"
#include "vars.h"

struct Parser::Imp {
  using Owner = Parser;

  Imp(Owner* owner, Vars& v) : owner_(owner), v_(v) {}

  static std::string RemoveComment(std::string s) {
    return s.substr(0, s.find('#', 0));
  }
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

  Owner* owner_;
  Vars& v_;
};

Parser::Parser(Vars& v) : imp(new Imp(this, v)) {}

Parser::~Parser() = default;

void Parser::Imp::Cmd(std::string s) {
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

void Parser::Imp::CmdSet(std::string s) {
  std::string cmd, type, key, val;
  std::stringstream b(s);
  b >> std::skipws;
  b >> cmd >> type >> key;
  char c;
  // Read first non-ws character
  b >> c;
  if (b.good()) {
    // Read remaining line
    std::getline(b, val);
    val = c + val;
  }

  v_.SetStr(type, key, val);
}

void Parser::Imp::CmdSetNext(std::string s) {
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

  if (!v_.Int.Contains(key)) {
    v_.Int.Set(key, 0);
  }
  int& id = v_.Int[key];
  v_.SetStr(type, key + std::to_string(id), val);
  ++id;
}

void Parser::Imp::CmdDel(std::string s) {
  std::string cmd, key;
  std::stringstream b(s);
  b >> cmd >> key;
  if (!v_.Del(key)) {
    throw std::runtime_error("CmdDel(): unknown variable '" + key + "'");
  }
}

void Parser::Imp::CmdInclude(std::string s) {
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

  owner_->RunAll(r);
}

void Parser::Run(std::string s) {
  imp->Cmd(s);
}

void Parser::RunNext(std::istream& in) {
  std::string s;
  std::getline(in, s);
  imp->Cmd(s);
}

void Parser::RunAll(std::istream& in) {
  while (in) {
    RunNext(in);
  }
}

template <class T>
void Parser::Print(const Vars::Map<T>& m, std::ostream& out) {
  for (auto it = m.cbegin(); it != m.cend(); ++it) {
    out << "set " << m.GetTypeName() << " " << it->first << " "
        << m.GetStr(it->first) << std::endl;
  }
}

void Parser::PrintAll(std::ostream& out) const {
  Print(imp->v_.String, out);
  Print(imp->v_.Int, out);
  Print(imp->v_.Double, out);
  Print(imp->v_.Vect, out);
}

void Parser::PrintAll() const {
  PrintAll(std::cout);
}
