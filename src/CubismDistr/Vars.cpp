#include <sstream>
#include <iostream>
#include <cassert>

#include "Vars.h"

template <>
Vars::Map<int>& Vars::Get<int>() {
  return Int;
}

template <>
Vars::Map<double>& Vars::Get<double>() {
  return Double;
}

template <>
Vars::Map<std::string>& Vars::Get<std::string>() {
  return String;
}

template <>
Vars::Map<std::vector<double>>& Vars::Get<std::vector<double>>() {
  return Vect;
}

template <>
std::string Vars::Map<std::string>::GetTypeName() const {
  return "string";
}

template <>
std::string Vars::Map<int>::GetTypeName() const {
  return "int";
}

template <>
std::string Vars::Map<double>::GetTypeName() const {
  return "double";
}

template <>
std::string Vars::Map<std::vector<double>>::GetTypeName() const {
  return "vect";
}


template <class T>
std::string Vars::Map<T>::Print(Key k) const {
  std::stringstream b;
  b << m_.at(k);
  return b.str();
}

template <class T>
bool Vars::Map<T>::Parse(std::string s, Key k) {
  std::stringstream b(s);
  b >> std::skipws >> m_[k];

  if (b.fail()) {
    std::cerr 
        << "Unable to parse '" << s 
        << "' as " << GetTypeName() << std::endl;;
    assert(false);
  }

  // Check that string contained only one element of type T
  // Read one more character
  char c;
  b >> std::skipws >> c;
  if (b.good()) {
    // If good, the string contained invalid characters
    std::cerr 
        << "Trailing characters when parsing '" << s 
        << "' as " << GetTypeName() << std::endl;;
    assert(false);
  }
  return true;
}

template <>
std::string Vars::Map<std::vector<double>>::Print(Key k) const {
  std::stringstream b;
  for (auto a : m_.at(k)) {
    b << a << " ";
  }
  return b.str();
}

template <>
bool Vars::Map<std::vector<double>>::Parse(std::string s, Key k) {
  std::vector<double> r;
  std::stringstream b(s);
  while (b) {
    double a;
    b >> std::skipws >> a;
    if (!b.fail()) {
      r.push_back(a);
    } else if (!b.eof()) {
      // String contained invalid characters
      std::cerr 
          << "Unable to parse '" << s 
          << "' as " << GetTypeName() << std::endl;;
      assert(false);
    }
  }
  m_[k] = r;
  return true;
}

template <>
bool Vars::Map<std::string>::Parse(std::string s, Key k) {
  m_[k] = s;
  return true;
}

template <class T>
const T& Vars::Map<T>::operator[](Key k) const {
  return m_.at(k);
}

template <class T>
T& Vars::Map<T>::operator[](Key k) {
  return m_.at(k);
}

template <class T>
const T* Vars::Map<T>::operator()(Key k) const {
  if (m_.count(k)) {
    return &m_.at(k);
  }
  return nullptr;
}

template <class T>
T* Vars::Map<T>::operator()(Key k) {
  if (m_.count(k)) {
    return &m_.at(k);
  }
  return nullptr;
}

template <class T>
void Vars::Map<T>::Set(Key k, const T& v) {
  m_[k] = v;
}

template <class T>
bool Vars::Map<T>::Exists(Key k) const {
  return m_.count(k);
}

template <class T>
void Vars::Map<T>::Del(Key k) {
  m_.erase(m_.find(k));
}

std::string Vars::Print(std::string type, Key k) const {
  if (type == String.GetTypeName()) {
    return String.Print(k);
  }
  if (type == Int.GetTypeName()) {
    return Int.Print(k);
  }
  if (type == Double.GetTypeName()) {
    return Double.Print(k);
  }
  if (type == Vect.GetTypeName()) {
    return Vect.Print(k);
  }
  std::cerr << "Vars::Print(): Unknown type=" << type << std::endl;
  assert(false);
  return "";
}

bool Vars::Parse(std::string s, std::string type, Key k) {
  if (type == String.GetTypeName()) {
    return String.Parse(s, k);
  }
  if (type == Int.GetTypeName()) {
    return Int.Parse(s, k);
  }
  if (type == Double.GetTypeName()) {
    return Double.Parse(s, k);
  }
  if (type == Vect.GetTypeName()) {
    return Vect.Parse(s, k);
  }
  std::cerr << "Vars::Parse(): Unknown type=" << type << std::endl;
  assert(false);
  return false;
}


namespace test_vars {

void TestSimple() {
  Vars v;
  v.Parse("asdf", "string", "a");
  v.Parse("5", "int", "b");
  v.Parse("5.5", "double", "c");
  v.Parse("1 2 3 4", "vect", "d");

  assert(v.String["a"] == "asdf");
  assert(v.Int["b"] == 5);
  assert(v.Double["c"] == 5.5);
  assert(v.Vect["d"] == std::vector<double>({1., 2., 3., 4.}));
  assert(v.Vect["d"] != std::vector<double>({2., 2., 3., 4.}));
}

template <class T>
void TestParse(std::string s) {
  Vars vars;
  std::string k = "key";
  auto& m = vars.Get<T>();
  m.Parse(s, k);
  auto p = m.Print(k);
  std::cerr
      << m.GetTypeName() << ": "
      << "'" << s << "' == '" << p << "'" << std::endl;
  assert(s == p);
}

void TestTypeName(std::string s, std::string type) {
  Vars vars;
  std::string k = "key";
  bool f = vars.Parse(s, type, k);
  assert(f);
  std::string p = vars.Print(type, k);
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
  m.Parse(s, k);
  assert(m.Exists(k));
  m.Del(k);
  assert(!m.Exists(k));
}

void Test() {
  std::cerr << "\ntest_vars::Test()" << std::endl;

  std::cerr << "TestSimple" << std::endl;
  TestSimple();

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

} // namespace test_vars

/*
template class Vars::Map<std::string>;
template class Vars::Map<int>;
template class Vars::Map<double>;
template class Vars::Map<std::vector<double>>;
*/
