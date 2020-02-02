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
  static bool IsNameChar(char c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') ||
           (c >= 'A' && c <= 'Z') || c == '_';
  }
  std::string GetStr(std::string key) const {
    auto type = v_.GetTypeName(key);
    if (type != "") {
      return v_.GetStr(type, key);
    }
    throw std::runtime_error("GetStr(): undefined variable '" + key + "'");
  }
  std::string ExpandVariables(std::string str) const {
    enum class S { normal, dollar, varname, expand, varname_braced};
    S s = S::normal; // state
    size_t i = 0;
    std::string res; // result
    std::string varname;
    try {
      while (i < str.size()) {
        char c = str[i];
        auto next = [&i]() { ++i; };
        switch (s) {
          case S::normal: {
            if (c == '$') {
              s = S::dollar;
              next();
            } else {
              res += c;
              next();
            }
            break;
          }
          case S::dollar: {
            if (c == '{') {
              s = S::varname_braced;
              next();
            } else if (IsNameChar(c)) {
              s = S::varname;
            } else {
              throw std::runtime_error(
                  std::string() + "S::dollar: invalid character '" + c + "'");
            }
            break;
          }
          case S::varname: {
            if (IsNameChar(c) && c != '_') {
              varname += c;
              next();
            } else {
              s = S::expand;
            }
            break;
          }
          case S::varname_braced: {
            if (IsNameChar(c)) {
              varname += c;
              next();
            } else if (c == '}') {
              s = S::expand;
              next();
            } else {
              throw std::runtime_error(
                  std::string() + "S::varname_braced: invalid character '" + c +
                  "'");
            }
            break;
          }
          case S::expand: {
            res += GetStr(varname);
            varname = "";
            s = S::normal;
            break;
          }
        }
      }
      if (s == S::varname && varname != "") {
        res += GetStr(varname);
      } else if (s == S::varname_braced) {
        throw std::runtime_error("unmatched brace {");
      }
    } catch (const std::runtime_error& e) {
      throw std::runtime_error(
          "ExpandVariables(): error while parsing '" + str + "'\n" + e.what());
    }
    return res;
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
  s = ExpandVariables(s);
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
