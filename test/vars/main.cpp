#include <sstream>
#include <iostream>
#include <cassert>

#include "Vars.h"

namespace simple {

void Simple() {
  Vars v;
  v.SetStr("string", "a", "asdf");
  v.SetStr("int", "b", "5");
  v.SetStr("double", "c", "5.5");
  v.SetStr("vect", "d", "1 2 3 4");

  assert(v.String["a"] == "asdf");
  assert(v.Int["b"] == 5);
  assert(v.Double["c"] == 5.5);
  assert(v.Vect["d"] == std::vector<double>({1., 2., 3., 4.}));
  assert(v.Vect["d"] != std::vector<double>({2., 2., 3., 4.}));
}

} // namespace simple

template <class T>
void TestParse(std::string s) {
  Vars vars;
  std::string k = "key";
  auto& m = vars.Get<T>();
  m.SetStr(k, s);
  auto p = m.GetStr(k);
  std::cerr
      << m.GetTypeName() << ": "
      << "'" << s << "' == '" << p << "'" << std::endl;
  assert(s == p);
}

void TestTypeName(std::string s, std::string type) {
  Vars vars;
  std::string k = "key";
  vars.SetStr(type, k, s);
  std::string p = vars.GetStr(type, k);
  std::cerr
      << type << ": "
      << "'" << s << "' == '" << p << "'" << std::endl;
  assert(s == p);
}

template <class T>
void TestPtr(const T& v /*value*/, const T& vo /*other*/) {
  Vars vars;
  std::string k = "key";
  auto& m = vars.Get<T>();
  m.Set(k, v);
  const auto& cm = m;
  assert(m[k] == v);
  assert(*m(k) == v);
  assert(cm[k] == v);
  assert(*cm(k) == v);

  m[k] = vo;
  assert(cm[k] == vo);

  *m(k) = v;
  assert(cm[k] == v);
}

template <class T>
void TestDel(std::string s) {
  Vars vars;
  std::string k = "key";
  auto& m = vars.Get<T>();
  m.SetStr(k, s);
  assert(m.Exists(k));
  m.Del(k);
  assert(!m.Exists(k));
}

void Test() {
  std::cerr << "\ntest_vars::Test()" << std::endl;

  std::cerr << "TestParse" << std::endl;
  TestParse<std::string>("asdf");
  TestParse<int>("123");
  TestParse<double>("123.456");
  TestParse<std::vector<double>>("1.2 3.4 5.6 ");

  std::cerr << "TestTypeName" << std::endl;
  TestTypeName("asdf", "string");
  TestTypeName("123", "int");
  TestTypeName("123.456", "double");
  TestTypeName("1.2 3.4 5.6 ", "vect");

  std::cerr << "TestPtr" << std::endl;
  TestPtr<std::string>("asdf", "qwer");
  TestPtr<int>(123, 456);
  TestPtr<double>(123.456, 456.123);
  TestPtr<std::vector<double>>({1.2, 3.4, 5.6}, {5.6, 3.4});

  std::cerr << "TestDel" << std::endl;
  TestDel<std::string>("asdf");
  TestDel<int>("123");
  TestDel<double>("123.456");
  TestDel<std::vector<double>>("1.2 3.4 5.6 ");

  std::cerr << "test_vars::Test() done\n" << std::endl;
}


int main() {
  simple::Simple();

  Test();
}
