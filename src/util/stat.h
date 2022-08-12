// Created by Petr Karnakov on 15.07.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <algorithm>
#include <functional>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "geom/mesh.h"
#include "solver/embed.h"
#include "util/logger.h"

// Selects a non-void type.
template <class...>
struct NonVoid;

template <class T>
struct NonVoid<T> {
  using type = T;
};

template <class T, class... TT>
struct NonVoid<T, TT...> {
  using type = typename std::conditional<
      std::is_same<T, void>::value, typename NonVoid<TT...>::type, T>::type;
};

// Deduces the return type of Func() called with Args, void if invalid.
// Case of no arguments.
template <class Func, class... Args>
struct ResultOf {
  template <class>
  static void Eval(...);
  template <class U>
  static auto Eval(decltype(std::declval<U>()())* r) {
    return *r;
  }
  using type = decltype(Eval<Func>(0));
};

// Case of one or more arguments.
template <class Func, class T, class... TT>
struct ResultOf<Func, T, TT...> {
  template <class...>
  static void Eval(...);
  template <class U, class... UU>
  static auto Eval(decltype(std::declval<Func>()( //
      std::declval<U>(), std::declval<UU>()...))* r) {
    return *r;
  }
  using type = decltype(Eval<T, TT...>(0));
};

// Deduces the return type of `Func::operator() const` if unambiguous.
template <class Func>
class ResultOfDeducedArgs {
 private:
  template <class...>
  static void Eval(...);
  template <class R, class... Args>
  static R Ret(R (Func::*)(Args...) const);
  template <class T>
  static auto Eval(decltype(Ret(&T::operator()))* r) {
    return *r;
  }

 public:
  using type = decltype(Eval<Func>(0));
};

template <class M>
class Stat {
 public:
  enum class Reduction { none, sum, min, max };
  enum class Type { scal, vect };
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  struct Entry {
    std::string name;
    std::string desc;
    Reduction op = Reduction::none;
    Type type = Type::scal;
    std::function<Scal()> func;
    std::function<Vect()> func_vect;
    bool enabled = true;
    bool hidden = false;
    bool derived =
        false; // computed from other entires, postponed to second pass
    Entry(
        std::string name_, std::string desc_, Reduction op_,
        std::function<Scal()> func_, bool hidden_, bool derived_)
        : name(name_)
        , desc(desc_)
        , op(op_)
        , type(Type::scal)
        , func(func_)
        , enabled(true)
        , hidden(hidden_)
        , derived(derived_) {}
    Entry(
        std::string name_, std::string desc_, Reduction op_,
        std::function<Vect()> func_vect_, bool hidden_, bool derived_)
        : name(name_)
        , desc(desc_)
        , op(op_)
        , type(Type::vect)
        , func_vect(func_vect_)
        , enabled(true)
        , hidden(hidden_)
        , derived(derived_) {}
  };

  Stat(M& m_, const Embed<M>* eb = nullptr) : m(m_), eb_(eb), vect(this) {}
  Stat(const Stat&) = delete;

  template <class Func>
  struct ResultOfMesh {
    using type = typename ResultOfDeducedArgs<Func>::type;
    /*
    using type = typename NonVoid<
        typename ResultOf<Func>::type, //
        typename ResultOf<Func, IdxCell, M>::type, //
        typename ResultOf<Func, IdxCell, Embed<M>>::type>::type;
        */
  };
  template <class T>
  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<T(IdxCell c, const M&)> func, bool hidden = false) {
    AddName(name);
    AddMeshLoop(name, desc, op, func, hidden, m);
  }
  template <class T>
  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<T(IdxCell c, const Embed<M>&)> func, bool hidden = false) {
    AddName(name);
    fassert(eb_, "Add: Can't add the entry to Stat created without Embed.");
    AddMeshLoop(name, desc, op, func, hidden, *eb_);
  }
  template <class T>
  void Add(
      std::string name, std::string desc, Reduction op, std::function<T()> func,
      bool hidden = false) {
    AddName(name);
    entries_.emplace(name, Entry(name, desc, op, func, hidden, false));
  }
  template <class Func>
  void AddDerived(
      std::string name, std::string desc, Func func, bool hidden = false) {
    AddName(name);
    entries_.emplace(
        name, Entry(
                  name, desc, Reduction::none,
                  [this, func]() { return func(*this); }, hidden, true));
  }
  template <class Func>
  void AddSum(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::sum, func);
  }
  template <class Func>
  void AddNone(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::none, func);
  }
  template <class Func>
  void AddMax(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::max, func);
  }
  template <class Func>
  void AddMin(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::min, func);
  }
  template <class Func>
  void AddSumHidden(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::sum, func, true);
  }
  template <class Func>
  void AddNoneHidden(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::none, func, true);
  }
  template <class Func>
  void AddMaxHidden(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::max, func, true);
  }
  template <class Func>
  void AddMinHidden(std::string name, std::string desc, Func func) {
    using Result = typename ResultOfMesh<Func>::type;
    Add<Result>(name, desc, Reduction::min, func, true);
  }
  void Update() {
    auto sem = m.GetSem();
    if (sem("reduce")) {
      for (auto& p : entries_) {
        auto& e = p.second;
        if (e.enabled && !e.derived) {
          switch (e.type) {
            case Type::scal:
              fassert(e.func);
              values_[e.name] = e.func();
              break;
            case Type::vect:
              fassert(e.func_vect);
              values_vect_[e.name] = e.func_vect();
              break;
          }
        }
      }
      for (auto& p : values_) {
        auto& e = entries_.at(p.first);
        if (e.enabled && !e.derived) {
          switch (e.op) {
            case Reduction::none:
              // nop
              break;
            case Reduction::sum:
              m.Reduce(&p.second, "sum");
              break;
            case Reduction::min:
              m.Reduce(&p.second, "min");
              break;
            case Reduction::max:
              m.Reduce(&p.second, "max");
              break;
          }
        }
      }
      for (auto& p : values_vect_) {
        auto& e = entries_.at(p.first);
        if (e.enabled && !e.derived) {
          for (size_t d = 0; d < dim; ++d) {
            switch (e.op) {
              case Reduction::none:
                // nop
                break;
              case Reduction::sum:
                m.Reduce(&p.second[d], "sum");
                break;
              case Reduction::min:
                m.Reduce(&p.second[d], "min");
                break;
              case Reduction::max:
                m.Reduce(&p.second[d], "max");
                break;
            }
          }
        }
      }
    }
    if (sem("derive")) {
      for (auto& p : entries_) {
        auto& e = p.second;
        if (e.enabled && e.derived) {
          switch (e.type) {
            case Type::scal:
              if (!values_.count(e.name)) {
                values_[e.name] = 0;
              }
              values_[e.name] = p.second.func();
              break;
            case Type::vect:
              if (!values_.count(e.name)) {
                values_vect_[e.name] = Vect(0);
              }
              values_vect_[e.name] = p.second.func_vect();
              break;
          }
        }
      }
    }
  }
  void WriteHeader(std::ostream& out, bool with_hidden = false) const {
    bool first = true;
    for (auto n : names_) {
      auto& e = entries_.at(n);
      if (e.enabled && (with_hidden || !e.hidden)) {
        first || out << ' ', first = false;
        switch (e.type) {
          case Type::scal:
            out << n;
            break;
          case Type::vect:
            for (size_t i : M::dirs) {
              if (i) out << ' ';
              out << n << '.' << M::direction(i).letter();
            }
            break;
        }
      }
    }
    out << std::endl;
  }
  void WriteValues(std::ostream& out, bool with_hidden = false) const {
    bool first = true;
    for (auto n : names_) {
      auto& e = entries_.at(n);
      if (e.enabled && (with_hidden || !e.hidden)) {
        first || out << ' ', first = false;
        switch (e.type) {
          case Type::scal:
            out << values_.at(n);
            break;
          case Type::vect:
            auto& v = values_vect_.at(n);
            for (size_t i : M::dirs) {
              if (i) out << ' ';
              out << v[i];
            }
            break;
        }
      }
    }
    out << std::endl;
  }
  void WriteSummary(std::ostream& out, bool with_hidden = false) const {
    for (auto n : names_) {
      if (with_hidden || !entries_.at(n).hidden) {
        out << n << ": " << entries_.at(n).desc << '\n';
      }
    }
  }
  void SortNames() {
    std::sort(names_.begin(), names_.end());
  }
  const std::vector<std::string>& GetNames() const {
    return names_;
  }
  bool Exists(std::string name) const {
    return entries_.count(name);
  }
  bool ExistsScal(std::string name) const {
    return entries_.count(name) && entries_[name].type == Type::scal;
  }
  bool ExistsVect(std::string name) const {
    return entries_.count(name) && entries_[name].type == Type::vect;
  }
  bool IsHidden(std::string name) const {
    return entries_[name].hidden;
  }
  bool IsEnabled(std::string name) const {
    return entries_[name].enabled;
  }
  bool SetEnabled(std::string name, bool enabled) {
    return entries_.at(name).enabled = enabled;
  }
  std::string GetDesc(std::string name) const {
    return entries_[name].desc;
  }
  Scal operator[](std::string name) const {
    fassert(values_.count(name), "entry '" + name + "' not found");
    return values_.at(name);
  }

 private:
  void CheckExisting(std::string name) const {
    fassert(!entries_.count(name), "Add: Entry '" + name + "' already exists.");
  }
  void AddName(std::string name) {
    CheckExisting(name);
    names_.push_back(name);
  }
  template <class T, class MEB>
  void AddMeshLoop(
      std::string name, std::string desc, Reduction op,
      std::function<T(IdxCell c, const MEB&)> func, bool hidden,
      const MEB& meb) {
    switch (op) {
      case Reduction::none:
        fassert(false, "Reduction::none not allowed here");
      case Reduction::sum:
        entries_.emplace(
            name, //
            Entry(
                name, desc, op,
                [func, &meb]() {
                  T sum(0);
                  for (auto c : meb.Cells()) {
                    sum += func(c, meb);
                  }
                  return sum;
                },
                hidden, false));
        break;
      case Reduction::max:
        entries_.emplace(
            name, //
            Entry(
                name, desc, op,
                [func, &meb]() {
                  T max(-std::numeric_limits<Scal>::max());
                  for (auto c : meb.Cells()) {
                    max = std::max(max, func(c, meb));
                  }
                  return max;
                },
                hidden, false));
        break;
      case Reduction::min:
        entries_.emplace(
            name, //
            Entry(
                name, desc, op,
                [func, &meb]() {
                  T min(std::numeric_limits<Scal>::max());
                  for (auto c : meb.Cells()) {
                    min = std::min(min, func(c, meb));
                  }
                  return min;
                },
                hidden, false));
        break;
    }
  }

  class LazyVect {
   public:
    LazyVect(const Stat* stat) : stat_(stat) {}
    Vect operator[](std::string name) const {
      fassert(
          stat_->values_vect_.count(name),
          "vect entry '" + name + "' not found");
      return stat_->values_vect_.at(name);
    }

   private:
    const Stat* const stat_;
  };

  M& m;
  const Embed<M>* eb_;
  std::vector<std::string> names_;
  std::map<std::string, Entry> entries_;
  std::map<std::string, Scal> values_;
  std::map<std::string, Vect> values_vect_;

 public:
  LazyVect vect;
};
