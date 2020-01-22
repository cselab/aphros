// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <map>
#include <string>
#include <vector>

class Vars {
 public:
  using Key = std::string;

  template <class T>
  class Map {
   public:
    using Value = T;
    using M = std::map<Key, Value>;
    using Iterator = typename M::iterator;
    using ConstIterator = typename M::const_iterator;
    std::string GetTypeName() const;
    std::string GetStr(Key) const;
    void SetStr(Key, std::string val);
    Value operator()(Key, Value) const;
    Value& operator[](Key);
    Value* Find(Key);
    const Value* Find(Key) const;
    const Value& operator[](Key) const;
    void Set(Key, const Value& val);
    bool Contains(Key) const;
    void Del(Key);
    bool DelIfContains(Key);
    Iterator begin() {
      return m_.begin();
    }
    Iterator end() {
      return m_.end();
    }
    ConstIterator cbegin() const {
      return m_.cbegin();
    }
    ConstIterator cend() const {
      return m_.cend();
    }

   private:
    M m_;
  };

  // Returns map by type
  template <class Value>
  Map<Value>& Get();
  template <class Value>
  const Map<Value>& Get() const;

  std::string GetStr(std::string type, Key) const;
  void SetStr(std::string type, Key, std::string val);
  // Returns type name of variable by key
  // (if multiple types found, returns one; if none, returns "")
  std::string GetTypeName(Key) const;
  // Deletes entry by key for all types, returns true if found
  bool Del(Key);

  Map<std::string> String;
  Map<int> Int;
  Map<double> Double;
  Map<std::vector<double>> Vect;
};
