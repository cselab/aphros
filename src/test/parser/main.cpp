// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#include <fstream>
#include <iostream>
#include <sstream>

#include "parse/argparse.h"
#include "parse/config.h"
#include "parse/parser.h"
#include "parse/vars.h"
#include "util/format.h"
#include "util/logger.h"

namespace simple {

void Simple() {
  std::cout << "\n" << __func__ << std::endl;
  Vars var;
  Parser parser(var);

  std::stringstream s;
  s << "set string a 1" << std::endl;
  s << "set int b 1" << std::endl;
  s << "set double c 1" << std::endl;
  s << "set vect d 1" << std::endl;

  parser.ParseStream(s);
  parser.PrintVars(std::cout);

  fassert_equal(var.String["a"], "1");
  fassert_equal(var.Int["b"], 1);
  fassert_equal(var.Double["c"], 1.);
  fassert(var.Vect["d"] == std::vector<double>({1}));

  // FindByKey
  fassert_equal(var.FindByKey("a").type, "string");
  fassert_equal(var.FindByKey("b").type, "int");
  fassert_equal(var.FindByKey("c").type, "double");
  fassert_equal(var.FindByKey("d").type, "vect");
  fassert(var.FindByKey("b").found);
  fassert_equal(var.FindByKey("b").value, "1");
  fassert_equal(var.FindByKey("b").key, "b");

  // GenerateSetCommand
  fassert_equal(Parser::GenerateSetCommand(var, "a", "aa"), "set string aa 1");
  fassert_equal(Parser::GenerateSetCommand(var, "a"), "set string a 1");
}

} // namespace simple

void TestFile() {
  std::cout << "\n" << __func__ << std::endl;
  Vars par;
  Parser parser(par);
  parser.ParseFile("a.conf");

  parser.PrintVars(std::cout);
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
  parser.AddVariable<int>({"--int", "-i"}, 3).Options({3, 4});
  parser.AddVariable<double>("--double", 3).Help("Type double").Options({3, 4});
  parser.AddVariable<std::string>("--string", "a");
  parser.AddVariable<std::vector<double>>("--vect", {0});
  parser.AddVariable<int>({"--int0", "-i0"});
  parser.AddVariable<double>("--double0");
  parser.AddVariable<std::string>("--string0");
  parser.AddVariable<std::vector<double>>({"--vect0", "-v0"})
      .Help("List of doubles");
  parser.AddVariable<int>("nx").Help("Size in x-direction");
  parser.AddVariable<int>("ny").Help("Size in y-direction");
  parser.AddVariable<int>("nz", 1).Help("Size in z-direction");
  parser.AddVariable<int>("bs", 32).Help("Block size").Options({8, 16, 32});
  parser.AddVariable<int>("bsy").Help("Block size in y").Options({8, 16, 32});

  std::cout << '\n';
  parser.PrintHelp(std::cout, true, "program");
  std::cout << '\n';

  std::cout << "\nKnown args:\n";
  parser.GetKnownArgs().ForEachMap([](const auto& map) {
    for (auto p : map) {
      std::cout << map.GetTypeName() << ' ' << p.first << ' '
                << map.GetStr(p.first) << '\n';
    }
  });

  std::cout << "\nParsed args:\n";
  auto args = parser.ParseArgs({
      "--double0", "5.4", "--double", "3.4", "--int", "4", "-i", "7",
      "16", // nx
      "8", // ny
      "8", // nz
      "9", // bs
      "8", // bsy
  });
  args.ForEachMap([](const auto& map) {
    for (auto p : map) {
      std::cout << map.GetTypeName() << ' ' << p.first << ' '
                << map.GetStr(p.first) << '\n';
    }
  });
}

void TestFormat() {
  std::cout << "\n" << __func__ << std::endl;

  int i = 5;
  std::string s = "abcdefgh";
  std::cout << util::Format(
      "Char '{0}' at position {1} in string '{2}'\n", s[i], i, s);
  std::cout << util::Format(
      "scientific {0:.3e}, fixed {0:.3f}, default {0:.3g}\n", M_PI);
  std::cout << util::Format(
      "scientific {:.3e}, fixed {:.3f}, default {:.3g}\n", M_PI, M_PI, M_PI);
  std::cout << util::Format(
      "scientific {:.3e}, fixed {:.3f}, default {:.3g}\n", 0.0, 1.1, 2.2);
  std::cout << util::Format(
      "scientific {0:.3e}, fixed {1:.3f}, default {2:.3g}\n", 0.0, 1.1, 2.2);
  std::cout << util::Format(
      "width scientific {0:10.3e}, fixed {0:10.3f}, default {0:10.3g}\n", M_PI);
  std::cout << util::Format(
      "width leadzero scientific {0:010.e}, fixed {0:010.3f}, default "
      "{0:010.3g}\n",
      M_PI);
  std::cout << util::Format(
      "leadzero scientific {0:010.e}, fixed {0:010.3f}, default {0:010.3g}\n",
      -M_PI);
  std::cout << util::Format(
      "width leadzero {:03d} {:05} {:08d}\n", 1, 123, 123456);
  using Vect = generic::Vect<double, 3>;
  std::cout << util::Format("Vect {} {0:.10f}\n", Vect(0., 1., 2.));
  std::cout << //
      util::Format(
          "a={0} b={0:} c={1:7} d={0:7.3f} e={0:07.3f} f={0:07.3g} g={0:.3e}\n",
          M_PI, 3.14);
}

int main() {
  simple::Simple();

  TestFile();
  TestConfig();
  TestArgumentParser();
  TestFormat();
}
