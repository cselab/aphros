#include <sstream>
#include <iostream>
#include <iostream>
#include <stdexcept>

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
std::string Vars::Map<T>::GetStr(Key k) const {
  std::stringstream b;
  b << m_.at(k);
  return b.str();
}

template <class T>
void Vars::Map<T>::SetStr(Key k, std::string v) {
  std::stringstream b(v);
  b >> std::skipws;
  b >> m_[k];

  if (b.fail()) {
    std::cerr 
        << "Unable to parse '" << v 
        << "' as " << GetTypeName() << std::endl;;
    throw std::runtime_error("Vars::SetStr(): Unable to parse");
  }

  // Check that string contained only one element of type T
  // Read one more character
  char c;
  b >> c;
  if (b.good()) {
    // If good, the string contained invalid characters
    std::cerr 
        << "Trailing characters when parsing '" << v
        << "' as " << GetTypeName() << std::endl;;
    throw std::runtime_error("Vars::GetStr(): Trailing characters");
  }
}

template <>
std::string Vars::Map<std::vector<double>>::GetStr(Key k) const {
  std::stringstream b;
  for (auto a : m_.at(k)) {
    b << a << " ";
  }
  return b.str();
}

template <>
void Vars::Map<std::vector<double>>::SetStr(Key k, std::string s) {
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
      throw std::runtime_error("Vars::GetStr(): Unable to parse");
    }
  }
  m_[k] = r;
}

template <>
void Vars::Map<std::string>::SetStr(Key k, std::string s) {
  m_[k] = s;
}

template <class T>
const T& Vars::Map<T>::operator[](Key k) const {
  if (!m_.count(k)) {
    std::cerr << "variable '" << k
        << "' of type '" << GetTypeName()
        << "' not found" << std::endl;
    throw std::runtime_error("Vars::GetStr(): Variable not found");
  }
  return m_.at(k);
}

template <class T>
T& Vars::Map<T>::operator[](Key k) {
  if (!m_.count(k)) {
    std::cerr << "variable '" << k
        << "' of type '" << GetTypeName()
        << "' not found" << std::endl;
    throw std::runtime_error("Vars::GetStr(): Variable not found");
  }
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

template <class T>
bool Vars::Map<T>::DelIfExists(Key k) {
  if (Exists(k)) {
    Del(k);
    return true;
  }
  return false;
}

bool Vars::Del(Key k) {
  bool r = false;
  r = r || String.DelIfExists(k);
  r = r || Int.DelIfExists(k);
  r = r || Double.DelIfExists(k);
  r = r || Vect.DelIfExists(k);
  return r;
}

std::string Vars::GetStr(std::string t, Key k) const {
  if (t == String.GetTypeName()) {
    return String.GetStr(k);
  }
  if (t == Int.GetTypeName()) {
    return Int.GetStr(k);
  }
  if (t == Double.GetTypeName()) {
    return Double.GetStr(k);
  }
  if (t == Vect.GetTypeName()) {
    return Vect.GetStr(k);
  }
  throw std::runtime_error("Vars:: GetStr(): Unknown type " + t);
  return "";
}

void Vars::SetStr(std::string t, Key k, std::string v) {
  std::string e = GetTypeName(k); // existing type
  if (e != "" && e != t) {
    std::cerr 
        << "Vars::SetStr(): Attempt to change type of variable '" 
        << k << "' from '" << t << "' to '" << e << "'" << std::endl;
    throw std::runtime_error("Vars::SetStr(): Attempt to change type");
  }
  // TODO: Map::Set() still allows same name for two variables

  if (t == String.GetTypeName()) {
    String.SetStr(k, v);
  } else if (t == Int.GetTypeName()) {
    Int.SetStr(k, v);
  } else if (t == Double.GetTypeName()) {
    Double.SetStr(k, v);
  } else if (t == Vect.GetTypeName()) {
    Vect.SetStr(k, v);
  } else {
    throw std::runtime_error("Vars:: SetStr(): Unknown type " + t);
  }
}

template class Vars::Map<std::string>;
template class Vars::Map<int>;
template class Vars::Map<double>;
template class Vars::Map<std::vector<double>>;
