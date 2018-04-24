#pragma once

#include <string>
#include <vector>
#include <map>

class Vars {
 public:
  using Key = std::string;

  template <class T>
  class Map {
   public:
    using ValueType = T;
    using M = std::map<Key, T>;
    using It = typename M::iterator;
    std::string GetTypeName() const;
    std::string GetStr(Key) const;
    void SetStr(Key, std::string val);
    T* operator()(Key);
    const T* operator()(Key) const;
    T& operator[](Key);
    const T& operator[](Key) const;
    void Set(Key, const T& val);
    bool Exists(Key) const;
    void Del(Key);
    bool DelIfExists(Key);
    It begin() { return m_.begin(); }
    It end() { return m_.end(); }

   private:
    M m_;
  };

  // Returns map by type
  template <class T>
  Map<T>& Get();

  std::string GetStr(std::string type, Key) const;
  void SetStr(std::string type, Key, std::string val);
  // Returns a type name for which entry with given k exists
  // (if multiple found, returns one; if none, returns "")
  std::string GetTypeName(Key) const;
  // Deletes entry by key for all types, returns true if found
  bool Del(Key); 

  Map<std::string> String;
  Map<int> Int;
  Map<double> Double;
  Map<std::vector<double>> Vect;
};

