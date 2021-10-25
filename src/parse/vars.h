// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <functional>
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

    Map() : hook_([](const Key&) {}) {}
    Map(std::function<void(const Key&)> hook) : hook_(hook) {}

    std::string GetTypeName() const;
    std::string GetStr(Key) const;
    void SetStr(Key, std::string val);
    Value operator()(Key, Value def) const;
    Value& operator[](Key);
    Value* Find(Key);
    const Value* Find(Key) const;
    const Value& operator[](Key) const;
    void Set(Key, const Value& val);
    bool Contains(Key) const;
    void Del(Key);
    bool DelIfContains(Key);
    auto begin() {
      return m_.begin();
    }
    auto end() {
      return m_.end();
    }
    auto begin() const {
      return m_.begin();
    }
    auto end() const {
      return m_.end();
    }
    auto cbegin() const {
      return m_.cbegin();
    }
    auto cend() const {
      return m_.cend();
    }
    static std::string ValueToStr(Value value);

   private:
    std::map<Key, Value> m_;
    std::function<void(const Key&)> hook_;
  };

  Vars() = default;
  Vars(std::function<void(const Key&)> hook)
      : String(hook), Int(hook), Double(hook), Vect(hook) {}

  // Returns map by type
  template <class Value>
  Map<Value>& Get();
  template <class Value>
  const Map<Value>& Get() const;

  std::string GetStr(std::string type, Key) const;
  void SetStr(std::string type, Key, std::string val);
  // Returns type name of variable by key.
  // If multiple types found, returns one; if none, returns "".
  std::string GetTypeName(Key) const;
  struct EntryAsString {
    bool found;
    std::string type;
    std::string key;
    std::string value;
  };
  EntryAsString FindByKey(std::string key) const;
  // Deletes entry by key for all types, returns true if found
  bool Del(Key);
  template <class F>
  void ForEachMap(F lambda) const {
    lambda(String);
    lambda(Int);
    lambda(Double);
    lambda(Vect);
  }

  Map<std::string> String;
  Map<int> Int;
  Map<double> Double;
  Map<std::vector<double>> Vect;
};
