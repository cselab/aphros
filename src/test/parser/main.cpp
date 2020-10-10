// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "parse/argparse.h"
#include "parse/config.h"
#include "parse/parser.h"
#include "parse/vars.h"
#include "util/logger.h"
#include "util/format.h"

namespace simple {

void Simple() {
  std::cout << "\n" << __func__ << std::endl;
  Vars par;
  Parser parser(par);

  std::stringstream s;
  s << "set string a 1" << std::endl;
  s << "set int b 1" << std::endl;
  s << "set double c 1" << std::endl;
  s << "set vect d 1" << std::endl;

  parser.ParseStream(s);
  parser.PrintAll(std::cout);

  assert(par.String["a"] == "1");
  assert(par.Int["b"] == 1);
  assert(par.Double["c"] == 1.);
  assert(par.Vect["d"] == std::vector<double>({1}));
}

} // namespace simple

void TestFile() {
  std::cout << "\n" << __func__ << std::endl;
  Vars par;
  Parser parser(par);
  parser.ParseFile("a.conf");

  // print variables to cout
  parser.PrintAll(std::cout);
}

struct Config : public ConfigBase {
  VAR_DOUBLE(height);
  VAR_INT(size);
  VAR_VECT(elems);
  VAR_VECT3(gravity);
  VAR_STRING(name);
  VAR_BOOL(enable_fluid);
};

void TestConfig() {
  std::cout << "\n" << __func__ << std::endl;

  using Vect = generic::Vect<double, 3>;

  Vars var;

  {
    Parser parser(var);
    std::stringstream s(R"EOF(
  set double height 1.2
  set int size 3
  set vect elems 1 2 3
  set vect gravity 4 5 6
  set string name name
  set int enable_fluid 3
  )EOF");
    parser.ParseStream(s);
  }

  Config config;
  config.Read(var);
  std::cout << config.height << std::endl;
  std::cout << config.size << std::endl;
  std::cout << Vect(config.elems) << std::endl;
  std::cout << config.gravity << std::endl;
  std::cout << config.name << std::endl;
  std::cout << config.enable_fluid << std::endl;
}

void TestArgumentParser() {
  std::cout << "\n" << __func__ << std::endl;

  ArgumentParser parser("Description");
  parser.AddVariable<int>({"--int", "-i"}, 3);
  parser.AddVariable<double>("--double", 3);
  parser.AddVariable<std::string>("--string", "a");
  parser.AddVariable<std::vector<double>>("--vect", {0});
  parser.AddVariable<int>({"--int0", "-i0"});
  parser.AddVariable<double>("--double0");
  parser.AddVariable<std::string>("--string0");
  parser.AddVariable<std::vector<double>>({"--vect0", "-v0"})
      .Help("List of doubles");
  parser.AddVariable<int>("nx").Help("size in x-direction");
  parser.AddVariable<int>("ny").Help("size in y-direction");
  parser.AddVariable<int>("nz", 1).Help("size in z-direction");

  std::cout << '\n';
  parser.PrintHelp(std::cout, true, "program");
  std::cout << '\n';

  std::cout << "\nKnown args:\n";
  parser.GetKnownArgs().ForEachMap([](const auto& map) {
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
      std::cout << map.GetTypeName() << ' ' << it->first << ' '
                << map.GetStr(it->first) << '\n';
    }
  });

  std::cout << "\nParsed args:\n";
  auto args = parser.ParseArgs({
      "--help", "--double0", "5.4", "-h", "--double", "2.4", "--int", "4", "-i",
      "7",
      "16", // nx
      "8", // ny
  });
  args.ForEachMap([](const auto& map) {
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
      std::cout << map.GetTypeName() << ' ' << it->first << ' '
                << map.GetStr(it->first) << '\n';
    }
  });
}

void TestFormat() {
  std::cout << "\n" << __func__ << std::endl;

  int i = 5;
  std::string s = "qwerasdf";
  std::cout << util::Format(
      "Char '{0}' at position {1} in string '{2}'\n", s[i], i, s);
  using Vect = generic::Vect<double, 3>;
  std::cout << util::Format("Vect {}\n", Vect(0., 1., 2.));
}


int main() {
  simple::Simple();

  TestFile();
  TestConfig();
  TestArgumentParser();
  TestFormat();
}
