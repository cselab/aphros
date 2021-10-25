// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "util/logger.h"
#include "vars.h"

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
const Vars::Map<int>& Vars::Get<int>() const {
  return Int;
}

template <>
const Vars::Map<double>& Vars::Get<double>() const {
  return Double;
}

template <>
const Vars::Map<std::string>& Vars::Get<std::string>() const {
  return String;
}

template <>
const Vars::Map<std::vector<double>>& Vars::Get<std::vector<double>>() const {
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
  if (String.Contains(k)) {
    return String.GetTypeName();
  }
  if (Int.Contains(k)) {
    return Int.GetTypeName();
  }
  if (Double.Contains(k)) {
    return Double.GetTypeName();
  }
  if (Vect.Contains(k)) {
    return Vect.GetTypeName();
  }
  return "";
}

template <class T>
std::string Vars::Map<T>::ValueToStr(Value value) {
  std::stringstream buf;
  buf << value;
  return buf.str();
}

template <>
std::string Vars::Map<std::vector<double>>::ValueToStr(Value value) {
  std::stringstream buf;
  bool first = true;
  for (auto a : value) {
    if (!first) {
      buf << ' ';
    } else {
      first = false;
    }
    buf << a;
  }
  return buf.str();
}

template <class T>
std::string Vars::Map<T>::GetStr(Key k) const {
  return ValueToStr(m_.at(k));
}

template <class T>
void Vars::Map<T>::SetStr(Key k, std::string v) {
  std::stringstream b(v);
  b >> std::skipws;
  b >> m_[k];

  fassert(
      !b.fail(), //
      "unable to parse '" + v + "' as '" + GetTypeName() +
          "' for variable named '" + k + "'");

  // Check that string contained only one element of type T
  // Read one more character
  char c;
  b >> c;
  fassert(
      !b.good(), //
      "trailing characters '" + v + "' as '" + GetTypeName() +
          "' for variable named '" + k + "'");
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
      fassert(
          false, //
          "unable to parse '" + s + "' as '" + GetTypeName() +
              "' for variable named '" + k + "'");
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
  auto it = m_.find(k);
  fassert(
      it != m_.end(),
      "variable '" + k + "' of type '" + GetTypeName() + "' not found");
  hook_(k);
  return it->second;
}

template <class T>
T& Vars::Map<T>::operator[](Key k) {
  auto it = m_.find(k);
  fassert(
      it != m_.end(), //
      "variable '" + k + "' of type '" + GetTypeName() + "' not found");
  hook_(k);
  return it->second;
}

template <class T>
auto Vars::Map<T>::Find(Key k) -> Value* {
  auto it = m_.find(k);
  if (it == m_.end()) {
    return nullptr;
  }
  hook_(k);
  return &it->second;
}

template <class T>
auto Vars::Map<T>::Find(Key k) const -> const Value* {
  auto it = m_.find(k);
  if (it == m_.end()) {
    return nullptr;
  }
  hook_(k);
  return &it->second;
}

template <class T>
auto Vars::Map<T>::operator()(Key k, Value def) const -> Value {
  auto it = m_.find(k);
  if (it == m_.end()) {
    return def;
  }
  hook_(k);
  return it->second;
}

template <class T>
void Vars::Map<T>::Set(Key k, const T& v) {
  m_[k] = v;
}

template <class T>
bool Vars::Map<T>::Contains(Key k) const {
  return m_.count(k);
}

template <class T>
void Vars::Map<T>::Del(Key k) {
  m_.erase(m_.find(k));
}

template <class T>
bool Vars::Map<T>::DelIfContains(Key k) {
  if (Contains(k)) {
    Del(k);
    return true;
  }
  return false;
}

bool Vars::Del(Key k) {
  bool r = false;
  r = r || String.DelIfContains(k);
  r = r || Int.DelIfContains(k);
  r = r || Double.DelIfContains(k);
  r = r || Vect.DelIfContains(k);
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
  fassert(false, "Vars:: GetStr(): Unknown type " + t);
  return "";
}

auto Vars::FindByKey(std::string key) const -> EntryAsString {
  EntryAsString res;
  res.found = false;
  res.key = key;
  auto check = [&](auto& map){
    if (map.Contains(key)) {
      res.found = true;
      res.type = map.GetTypeName();
      res.value = map.GetStr(key);
    }
  };
  check(String);
  check(Int);
  check(Double);
  check(Vect);
  return res;
}

void Vars::SetStr(std::string t, Key k, std::string v) {
  std::string e = GetTypeName(k); // existing type
  if (e != "" && e != t) {
    fassert(
        false, //
        "attempt to change type of variable '" + k + "' from '" + e + "' to '" +
            t + "'");
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
    fassert(false, "unknown type " + t);
  }
}

template class Vars::Map<std::string>;
template class Vars::Map<int>;
template class Vars::Map<double>;
template class Vars::Map<std::vector<double>>;
