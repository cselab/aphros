#pragma once

#include <algorithm>
#include <functional>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "geom/mesh.h"
#include "solver/embed.h"
#include "util/logger.h"

template <class M>
class Stat {
 public:
  enum class Reduction { none, sum, min, max };
  using Scal = typename M::Scal;

  struct Entry {
    std::string name;
    std::string desc;
    Reduction op;
    std::function<Scal()> func;
    bool hidden;
    bool derived; // computed from other entires, postponed to second pass
  };

  Stat(M& m, const Embed<M>* eb = nullptr) : m(m), eb_(eb) {}
  Stat(const Stat&) = delete;

  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<Scal(IdxCell c, const M&)> func, bool hidden=false) {
    AddName(name);
    AddMeshLoop(name, desc, op, func, hidden, m);
  };
  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<Scal(IdxCell c, const Embed<M>&)> func, bool hidden=false) {
    AddName(name);
    fassert(eb_, "Add: Can't add the entry to Stat created without Embed.");
    AddMeshLoop(name, desc, op, func, hidden, *eb_);
  };
  void Add(
      std::string name, std::string desc, Reduction op,
      std::function<Scal()> func, bool hidden=false) {
    AddName(name);
    entries_.emplace(name, Entry{name, desc, op, func, hidden, false});
  };
  void AddDerived(
      std::string name, std::string desc, std::function<Scal(const Stat&)> func,
      bool hidden = false) {
    AddName(name);
    entries_.emplace(
        name, Entry{
                  name, desc, Reduction::none,
                  [this, func]() { return func(*this); }, hidden, true});
  };
  template <class Func>
  void AddSum(std::string name, std::string desc, Func func) {
    Add(name, desc, Reduction::sum, func);
  }
  template <class Func>
  void AddMax(std::string name, std::string desc, Func func) {
    Add(name, desc, Reduction::max, func);
  }
  template <class Func>
  void AddMin(std::string name, std::string desc, Func func) {
    Add(name, desc, Reduction::min, func, true);
  }
  template <class Func>
  void AddSumHidden(std::string name, std::string desc, Func func) {
    Add(name, desc, Reduction::sum, func, true);
  }
  template <class Func>
  void AddMaxHidden(std::string name, std::string desc, Func func) {
    Add(name, desc, Reduction::max, func, true);
  }
  template <class Func>
  void AddMinHidden(std::string name, std::string desc, Func func) {
    Add(name, desc, Reduction::min, func, true);
  }
  void Update() {
    auto sem = m.GetSem();
    if (sem("reduce")) {
      for (auto& p : entries_) {
        if (!p.second.derived) {
          values_[p.first] = p.second.func();
        }
      }
      for (auto& p : values_) {
        auto& e = entries_.at(p.first);
        if (!e.derived) {
          switch (e.op) {
            case Reduction::none:
              throw std::runtime_error(
                  FILELINE + ": Reduction::none not allowed here");
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
    }
    if (sem("derive")) {
      for (auto& p : entries_) {
        if (p.second.derived) {
          values_[p.first] = p.second.func();
        }
      }
    }
  }
  void WriteHeader(std::ostream& out, bool with_hidden=false) const {
    bool first = true;
    for (auto n : names_) {
      if (with_hidden || !entries_.at(n).hidden) {
        first || out << ' ', first = false;
        out << n;
      }
    }
    out << std::endl;
  }
  void WriteValues(std::ostream& out, bool with_hidden=false) const {
    bool first = true;
    for (auto n : names_) {
      if (with_hidden || !entries_.at(n).hidden) {
        first || out << ' ', first = false;
        out << values_.at(n);
      }
    }
    out << std::endl;
  }
  void WriteSummary(std::ostream& out, bool with_hidden=false) const {
    for (auto n : names_) {
      if (with_hidden || !entries_.at(n).hidden) {
        out << n << ": " << entries_.at(n).desc << '\n';
      }
    }
  }
  void SortNames() {
    std::sort(names_.begin(), names_.end());
  }
  const Entry& GetEntry(std::string name) const {
    return entries_[name];
  }
  const std::vector<std::string>& GetNames() const {
    return names_;
  }
  bool IsHidden(std::string name) const {
    return entries_[name].hidden;
  }
  std::string GetDesc(std::string name) const {
    return entries_[name].desc;
  }
  Scal operator[](std::string name) const {
    auto it = values_.find(name);
    if (it == values_.end()) {
      throw std::runtime_error(FILELINE + ": entry '" + name + "' not found");
    }
    return it->second;
  }

 private:
  void CheckExisting(std::string name) const {
    fassert(!entries_.count(name), "Add: Entry '" + name + "' already exists.");
  }
  void AddName(std::string name) {
    CheckExisting(name);
    names_.push_back(name);
  }
  template <class MEB>
  void AddMeshLoop(
      std::string name, std::string desc, Reduction op,
      std::function<Scal(IdxCell c, const MEB&)> func, bool hidden,
      const MEB& meb) {
    switch (op) {
      case Reduction::none:
        throw std::runtime_error(
            FILELINE + ": Reduction::none not allowed here");
      case Reduction::sum:
        entries_.emplace(
            name, //
            Entry{name, desc, op, [func, &meb = meb]() {
                    Scal sum = 0;
                    for (auto c : meb.Cells()) {
                      sum += func(c, meb);
                    }
                    return sum;
                  }, hidden, false});
        break;
      case Reduction::max:
        entries_.emplace(
            name, //
            Entry{name, desc, op, [func, &meb = meb]() {
                    Scal max = -std::numeric_limits<Scal>::max();
                    for (auto c : meb.Cells()) {
                      max = std::max(max, func(c, meb));
                    }
                    return max;
                  }, hidden, false});
        break;
      case Reduction::min:
        entries_.emplace(
            name, //
            Entry{name, desc, op, [func, &meb = meb]() {
                    Scal min = std::numeric_limits<Scal>::max();
                    for (auto c : meb.Cells()) {
                      min = std::min(min, func(c, meb));
                    }
                    return min;
                  }, hidden, false});
        break;
    }
  }

  M& m;
  const Embed<M>* eb_;
  std::vector<std::string> names_;
  std::map<std::string, Entry> entries_;
  std::map<std::string, Scal> values_;
};
