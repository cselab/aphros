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

std::string Vars::GetTypeName(Key k) const {
  if (String.Exists(k)) {
    return String.GetTypeName();
  }
  if (Int.Exists(k)) {
    return Int.GetTypeName();
  }
  if (Double.Exists(k)) {
    return Double.GetTypeName();
  }
  if (Vect.Exists(k)) {
    return Vect.GetTypeName();
  }
  return "";
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
  b >> std::skipws;
  b >> m_[k];

  if (b.fail()) {
    std::cerr 
        << "Unable to parse '" << s 
        << "' as " << GetTypeName() << std::endl;;
    assert(false);
  }

  // Check that string contained only one element of type T
  // Read one more character
  char c;
  b >> c;
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
  std::string e = GetTypeName(k); // existing type
  if (e != "" && e != type) {
    std::cerr << "Vars::Parse(): Attempt to change type of variable '" 
        << k << "' from '" << type << "' to '" << e << "'" << std::endl;
    assert(false);
    return false;
  }
  // TODO: Map::Set() still allows same name for two variables

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

template class Vars::Map<std::string>;
template class Vars::Map<int>;
template class Vars::Map<double>;
template class Vars::Map<std::vector<double>>;
