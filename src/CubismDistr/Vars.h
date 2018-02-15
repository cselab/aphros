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
    std::string Print(Key k) const;
    bool Parse(std::string s, Key k);
    T* operator()(Key k);
    const T* operator()(Key k) const;
    T& operator[](Key k);
    const T& operator[](Key k) const;
    void Set(Key k, const T& v);
    bool Exists(Key k) const;
    void Del(Key k);
    It begin() { return m_.begin(); }
    It end() { return m_.end(); }

   private:
    M m_;
  };

  // Returns map by type
  template <class T>
  Map<T>& Get();

  // Returns string representation by type name and key
  std::string Print(std::string type, Key k) const;
  // Parses string and sets value of given type and key
  bool Parse(std::string s, std::string type, Key k);

  Map<std::string> String;
  Map<int> Int;
  Map<double> Double;
  Map<std::vector<double>> Vect;
};

namespace test_vars {

void Test();

} // namespace test_vars
