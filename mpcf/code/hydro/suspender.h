#pragma once

#include <list>
#include <string>
#include <cassert>
#include <sstream>

// Suspendable functions.
// Function F() is separated in stages each enclosed by if-operator.
// At each call of function F(), only one stage is executed.
// Functions with stages can call other functions with stages
// in a separate stage.
class Suspender {
 public:
  struct U { // stage co[u]nter
    int c; // current
    int t; // target
    U(int c, int t) : c(c), t(t) {}
  };
  class Sem { // [sem]aphore
   public:
    // Constructor
    // Advance list iterator, add new counter if needed, reset counter
    Sem(Suspender& p, std::string name="") 
    : p(p), name_(name)
    {
      auto& l = p.lu_;
      auto& i = p.lui_;
      if (std::next(i) == l.end()) {
        l.emplace_back(0, 0);
      }
      ++i;
      i->c = 0;
    }
    // Forbid copy
    Sem(Sem&) = delete;
    Sem& operator=(Sem&) = delete;
    // Allow move
    Sem(Sem&&) = default;
    // Destructor
    // If all lower levels done, next stage.
    // If all stages on current level done, remove current level
    ~Sem() {
      auto& l = p.lu_;
      auto& i = p.lui_;

      assert(!l.empty());
      assert(i != l.end());
      assert(i != l.begin());

      auto ip = std::prev(i);

      if (std::next(i) == l.end()) {
        // all lower levels done, next stage
        ++i->t;
        if (i->c == i->t) { 
          // all stages done, remove current level
          // i->c keeps number of stages
          l.pop_back();
        }
      } 
      i = ip;
    }
    // Returns true if current stage needs execution
    // and advances stage counter
    bool operator()() {
      auto& i = p.lui_;
      return i->c++ == i->t;
    }
    std::string GetName() const {
      return name_;
    }
    void SetName(std::string name) {
      name_ = name;
    }
   private:
    Suspender& p; // parent
    std::string name_;
  };
  friend Sem;
  // Intializes list with auxiliary counter (-1,-1), sets iterator to it
  Suspender() 
    : lu_(1, U(-1,-1)), lui_(lu_.begin()) 
  {}
  Sem GetSem(std::string name="") {
    return Sem(*this, name);
  }
  // Converts counter list to string
  std::string LeToStr() const {
    std::stringstream b;
    for (auto e : lu_) {
      b << "(" << e.c << " " << e.t << ") ";
    }
    return b.str();
  }
  // Returns true if there are unfinished levels 
  bool Pending() const {
    return lu_.size() != 1;
  }

 private:
  using LU = std::list<U>;
  LU lu_;      // [l]ist of co[u]nters
  LU::iterator lui_; // [l]ist of co[u]nters [i]terator
};
