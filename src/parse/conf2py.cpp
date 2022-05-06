// Created by Petr Karnakov on 17.12.2020
// Copyright 2020 ETH Zurich

#include <fstream>
#include <iostream>
#include <memory>

#include "parse/argparse.h"
#include "parse/parser.h"
#include "parse/vars.h"

int main(int argc, const char** argv) {
  ArgumentParser argparser("Converts configuration to Python");
  argparser.AddVariable<std::string>("input", "-")
      .Help(
          "Path to configuration file with commands 'set ...'. If set to '-', "
          "read from stdin");
  argparser.AddVariable<std::string>("output", "-")
      .Help("Path to python output. If set to '-', write to stdout");

  auto args = argparser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  Vars var;
  Parser parser(var);
  const auto ipath = args.String["input"];
  if (ipath == "-") {
    parser.ParseStream(std::cin);
  } else {
    parser.ParseFile(ipath);
  }

  auto float_to_str = [](double a) {
    std::string s = util::Format("{:.16g}", a);
    if (s.find_first_of('.') == std::string::npos &&
        s.find_first_of('e') == std::string::npos) {
      s += '.';
    }
    return s;
  };

  const auto opath = args.String["output"];
  std::ofstream fout;
  if (opath != "-") {
    fout.open(opath);
  }
  auto& out = (opath != "-" ? fout : std::cout);
  for (auto a : var.String) {
    if (a.second.find('\n') != std::string::npos) {
      out << util::Format("{} = \"\"\"{}\"\"\"\n", a.first, a.second);
    } else {
      out << util::Format("{} = \"{}\"\n", a.first, a.second);
    }
  }
  for (auto a : var.Int) {
    out << util::Format("{} = {:}\n", a.first, a.second);
  }
  for (auto a : var.Double) {
    out << util::Format("{} = {}\n", a.first, float_to_str(a.second));
  }
  for (auto a : var.Vect) {
    std::string val;
    bool first = true;
    for (auto q : a.second) {
      auto sep = first ? "" : ", ";
      val += util::Format("{}{}", sep, float_to_str(q));
      first = false;
    }
    out << util::Format("{} = [{}]\n", a.first, val);
  }
}
