// Created by Petr Karnakov on 07.08.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "evalexpr.h"
#include "parser.h"
#include "util/filesystem.h"
#include "util/format.h"
#include "util/logger.h"
#include "vars.h"

static std::string Strip(std::string s) {
  size_t b = 0;
  size_t e = s.size();
  while (b < e && std::isspace(s[b])) {
    ++b;
  }
  while (b < e && std::isspace(s[e - 1])) {
    --e;
  }
  return s.substr(b, e - b);
}

struct Parser::Imp {
  using Owner = Parser;

  Imp(Owner* owner, Vars& var) : owner_(owner), var_(var) {}

  static std::string RemoveComment(std::string s) {
    return s.substr(0, s.find('#', 0));
  }
  static bool IsNameChar(char c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') ||
           (c >= 'A' && c <= 'Z') || c == '_';
  }
  std::string GetStr(std::string key) const {
    auto type = var_.GetTypeName(key);
    if (type != "") {
      return var_.GetStr(type, key);
    }
    fassert(false, "undefined variable '" + key + "'");
  }
  std::string ExpandVariables(std::string str) const {
    enum class S {
      normal,
      dollar,
      escape,
      varname,
      expand,
      varname_braced,
      stop
    };
    S s = S::normal; // State.
    size_t i = 0;
    std::string res; // Result.
    std::string varname;
    const char kEnd = 0;
    try {
      while (s != S::stop) {
        char c = (i < str.size() ? str[i] : kEnd);
        auto next = [&i]() { ++i; };
        switch (s) {
          case S::normal: {
            if (c == '\\') {
              s = S::escape;
              next();
            } else if (c == '$') {
              s = S::dollar;
              next();
            } else if (c == kEnd) {
              s = S::stop;
            } else {
              res += c;
              next();
            }
            break;
          }
          case S::escape: {
            res += '\\';
            res += c;
            next();
            s = S::normal;
            break;
          }
          case S::dollar: {
            if (c == '{') {
              s = S::varname_braced;
              next();
            } else if (c == '(') {
              res += '$';
              s = S::normal;
            } else if (IsNameChar(c)) {
              s = S::varname;
            } else {
              throw std::runtime_error(util::Format(
                  "unexpected character '{}' following $ while parsing '{}'", c,
                  str));
            }
            break;
          }
          case S::varname: {
            if (IsNameChar(c)) {
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
            } else if (c == kEnd) {
              throw std::runtime_error(
                  util::Format("unmatched brace {{ while parsing '{}'", str));
            } else {
              throw std::runtime_error(util::Format(
                  "unexpected character '{}' inside {{}} while parsing '{}'", c,
                  str));
            }
            break;
          }
          case S::expand: {
            res += GetStr(varname);
            varname.clear();
            s = S::normal;
            break;
          }
          case S::stop:
            break;
        }
      }
    } catch (const std::runtime_error& e) {
      fassert(false, "error while parsing '" + str + "'\n" + e.what());
    }
    return res;
  }
  std::string EvalExpressions(std::string str) const {
    // Evaluates arithmetic expressions inside $((...)).
    enum class S {
      normal,
      dollar,
      escape,
      paren1,
      expr,
      parenclose1,
      evalexpr,
      stop
    };
    S s = S::normal; // State.
    size_t i = 0;
    int depth = 0; // Depth of nested parentheses inside $((...)).
    std::string res; // Result.
    std::string expr;
    const char kEnd = 0;
    try {
      while (s != S::stop) {
        char c = (i < str.size() ? str[i] : kEnd);
        auto next = [&i]() { ++i; };
        switch (s) {
          case S::normal: {
            if (c == '\\') {
              s = S::escape;
              next();
            } else if (c == '$') {
              s = S::dollar;
              next();
            } else if (c == kEnd) {
              s = S::stop;
            } else {
              res += c;
              next();
            }
            break;
          }
          case S::escape: {
            res += c;
            next();
            s = S::normal;
            break;
          }
          case S::dollar: {
            if (c == '(') {
              s = S::paren1;
              next();
            } else {
              res += '$';
              res += c;
              s = S::normal;
              next();
            }
            break;
          }
          case S::paren1: {
            if (c == '(') {
              s = S::expr;
              depth = 0;
              expr.clear();
              next();
            } else {
              throw std::runtime_error(util::Format(
                  "unexpected character '{}' following ( while parsing '{}'", c,
                  str));
            }
            break;
          }
          case S::expr: {
            if (c == ')') {
              if (depth > 0) {
                --depth;
                expr += c;
              } else {
                s = S::parenclose1;
              }
              next();
            } else if (c == '(') {
              ++depth;
              expr += c;
              next();
            } else if (c == kEnd) {
              throw std::runtime_error(util::Format(
                  "unexpected end of string while parsing '{}'", str));
            } else {
              expr += c;
              next();
            }
            break;
          }
          case S::parenclose1: {
            if (c == ')') {
              s = S::evalexpr;
              next();
            } else {
              throw std::runtime_error(util::Format(
                  "unexpected character '{}' after ')' inside $(()) while "
                  "parsing '{}'",
                  c, str));
            }
            break;
          }
          case S::evalexpr: {
            res += util::Format("{:.14g}", EvalExpr<double>(expr));
            s = S::normal;
            break;
          }
          case S::stop:
            break;
        }
      }
    } catch (const std::runtime_error& e) {
      fassert(false, "error while parsing '" + str + "'\n" + e.what());
    }
    return res;
  }
  // Executes single command
  // In cases of error, if `filename != ""` and `line >= 0`
  // appends the error message with curpath and line.
  void Cmd(std::string, std::string curpath, int line);
  // set <type> <key> <value>
  void CmdSet(std::string);
  // setnext <type> <key> <value>
  // Sets value with name key+id and increments id stored in var_.Int[key]
  void CmdSetNext(std::string);
  // del <name>
  void CmdDel(std::string);
  // include <path>
  // curpath: path to file from which the command is called
  void CmdInclude(std::string path, std::string curpath);

  Owner* owner_;
  Vars& var_;
};

Parser::Parser(Vars& var) : imp(new Imp(this, var)) {}

Parser::~Parser() = default;

void Parser::Imp::Cmd(std::string s, std::string curpath, int line) {
  try {
    s = RemoveComment(s);
    s = ExpandVariables(s);
    s = EvalExpressions(s);
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
      CmdInclude(s, curpath);
    } else if (cmd == "eval") {
      std::stringstream rest;
      rest << b.rdbuf();
      Cmd(rest.str(), curpath, line);
    } else if (cmd == "") {
      // nop
    } else {
      throw std::runtime_error("unknown command '" + cmd + "'");
    }
  } catch (const std::runtime_error& e) {
    std::string msg = e.what();
    if (curpath != "" && line >= 0) {
      const auto loc = util::GetRealpath(curpath) + ':' + std::to_string(line);
      msg = loc + ": required from here\n" + msg;
    }
    throw std::runtime_error(msg);
  }
}

void Parser::Imp::CmdSet(std::string s) {
  std::string cmd, type, key, val;
  std::stringstream b(s);
  b >> std::skipws;
  b >> cmd >> type >> key;
  b >> std::noskipws;
  char c = ' ';
  while (b && std::isspace(c)) {
    b >> c;
  }
  while (b) {
    val += c;
    b >> c;
  }
  val = Strip(val);
  var_.SetStr(type, key, val);
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

  if (!var_.Int.Contains(key)) {
    var_.Int.Set(key, 0);
  }
  int& id = var_.Int[key];
  var_.SetStr(type, key + std::to_string(id), val);
  ++id;
}

void Parser::Imp::CmdDel(std::string s) {
  std::string cmd, key;
  std::stringstream b(s);
  b >> cmd >> key;
  fassert(var_.Del(key), "CmdDel(): unknown variable '" + key + "'");
}

void Parser::Imp::CmdInclude(std::string s, std::string curpath) {
  std::string cmd;
  std::string filename;
  std::stringstream b(s);
  b >> cmd >> filename;
  std::string dir;
  owner_->ParseFile(filename, curpath != "" ? util::GetDirname(curpath) : "");
}

void Parser::Run(std::string s) {
  imp->Cmd(s, "", -1);
}

// Reads multiline string from stream.
// String is read until EOL or while inside quotation marks "".
// Removes quotation marks and strips whitespaces.
// Returns:
//   string read,
//   number of lines read from stream.
std::pair<std::string, int> ReadMultiline(std::istream& in) {
  int lines = 0;
  const auto flags = in.flags();
  in >> std::noskipws;
  std::string res;
  enum class S { normal, quote, exit };
  S s = S::normal; // state
  char c = ' ';
  try {
    while (in && s != S::exit) {
      auto next = [&in, &c, &lines]() {
        in >> c;
        if (c == '\n') {
          ++lines;
        }
      };
      switch (s) {
        case S::normal: {
          if (c == '"') {
            s = S::quote;
            next();
          } else if (c == '\n') {
            s = S::exit;
          } else {
            res += c;
            next();
          }
          break;
        }
        case S::quote: {
          if (c == '"') {
            s = S::normal;
          } else {
            res += c;
          }
          next();
          break;
        }
        case S::exit: {
          break;
        }
      }
    }
    if (s == S::quote) {
      throw std::runtime_error("no matching '\"'");
    }
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(
        std::string(__func__) + ": " + e.what() + "\nstopping at string: '" +
        Strip(res) + "'");
  }
  in.flags(flags);
  return std::make_pair(Strip(res), lines);
}

void Parser::RunNext(std::istream& in) {
  const std::pair<std::string, int> p = ReadMultiline(in);
  imp->Cmd(p.first, "", p.second);
}

void Parser::ParseStream(std::istream& in) {
  while (in) {
    RunNext(in);
  }
}

void Parser::ParseFile(std::string path, std::string dir) {
  fassert(path != "", "Empty path");
  std::ifstream f(path);
  if (dir != "") {
    fassert(
        util::IsDir(dir),
        util::Format(
            "Not a directory '{}' specified to look for file '{}'", dir, path));
    auto path2 = util::Join(dir, path);
    if (!f.good() && path[0] != '/') {
      f.open(path2);
      fassert(
          f.good(), util::Format("Can't open file '{}' or '{}'", path, path2));
    }
    path = path2;
  }
  fassert(f.good(), util::Format("Can't open file '{}'", path));

  int line = 0;
  while (f) {
    const std::pair<std::string, int> p = ReadMultiline(f);
    line += p.second;
    imp->Cmd(p.first, path, line);
  }
}

template <class T>
void Parser::PrintMap(const Vars::Map<T>& m, std::ostream& out) {
  for (auto it = m.cbegin(); it != m.cend(); ++it) {
    const std::string val = m.GetStr(it->first);
    out << "set " << m.GetTypeName() << " " << it->first << " ";
    if (val.find("\n") != std::string::npos) {
      out << "\"" << val << "\n\"";
    } else {
      out << val;
    }
    out << std::endl;
  }
}

void Parser::PrintVars(const Vars& var, std::ostream& out) {
  PrintMap(var.String, out);
  PrintMap(var.Int, out);
  PrintMap(var.Double, out);
  PrintMap(var.Vect, out);
}

void Parser::PrintVars(std::ostream& out) const {
  PrintVars(imp->var_, out);
}
