#include <string>
#include <vector>
#include <map>
#include <sstream>

class Vars {
 public:
  template <class T>
  class Map {
   public:
    using Key = std::string;
    std::string Print(Key k) const;
    bool Parse(std::string s, Key k);
    T* operator()(Key k);
    const T* operator()(Key k) const;
    T& operator[](Key k);
    const T& operator[](Key k) const;
    void Set(Key k, const T& v);
    bool Exists(Key k);
   private:
    std::map<Key, T> m_;
  };

  // Returns map by type
  template <class T>
  Map<T>& Get();

  Map<std::string> Str;
  Map<int> Int;
  Map<double> Double;
  Map<std::vector<double>> Vect;
};

template <class T>
std::string Vars::Map<T>::Print(Key k) const {
  std::stringstream b;
  b << m_.at(k);
  return b.str();
}

template <class T>
bool Vars::Map<T>::Parse(std::string s, Key k) {
  std::stringstream b(s);
  b >> m_[k];
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
    b >> a;
    if (!b.fail()) {
      r.push_back(a);
    }
  }
  m_[k] = r;
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
bool Vars::Map<T>::Exists(Key k) {
  return m_.count(k);
}

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
  return Str;
}

template <>
Vars::Map<std::vector<double>>& Vars::Get<std::vector<double>>() {
  return Vect;
}

namespace test_vars {

template <class T>
void TestParse(std::string s, std::string title="", bool fatal=true) {
  Vars vars;
  std::string k = "key";
  auto& m = vars.Get<T>();
  m.Parse(s, k);
  auto p = m.Print(k);
  std::cerr
      << title << ": "
      << "'" << s << "' == '" << p << "'" << std::endl;
  if (fatal) {
    assert(s == p);
  }
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

void Test() {
  TestParse<std::string>("asdf", "str");
  TestParse<int>("123", "int");
  TestParse<double>("123.456", "double");
  TestParse<std::vector<double>>("1.2 3.4 5.6 ", "vector<double>");

  TestPtr<std::string>("asdf", "qwer");
  TestPtr<int>(123, 456);
  TestPtr<double>(123.456, 456.123);
  TestPtr<std::vector<double>>({1.2, 3.4, 5.6}, {5.6, 3.4});
}

} // namespace test_vars
