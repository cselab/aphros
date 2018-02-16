#pragma once

#include <list>
#include <string>
#include <sstream>

// TODO: Sequence like 
//   GetSem();
//   sem = GetSem();
//   sem();
// causes an infinite loop
// (see todo in tests/suspender)

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
    Sem(Suspender& p, std::string name="");
    // Forbid copy
    Sem(Sem&) = delete;
    Sem& operator=(Sem&) = delete;
    // Allow move
    Sem(Sem&&) = default;
    // Destructor
    // If all lower levels done, next stage.
    // If all stages on current level done, remove current level
    ~Sem();
    // Returns true if current stage needs execution
    // and advances stage counter
    bool operator()();
    std::string GetName() const { return name_; }
    void SetName(std::string name) { name_ = name; }
   private:
    Suspender& p; // parent
    std::string name_;
  };
  friend Sem;
  // Intializes list with auxiliary counter (-1,-1), sets iterator to it
  Suspender();
  Sem GetSem(std::string name="");
  // Converts counter list to string
  std::string Print() const;
  // Returns true if there are unfinished levels 
  bool Pending() const;

 private:
  using LU = std::list<U>;
  LU lu_;      // [l]ist of co[u]nters
  LU::iterator lui_; // [l]ist of co[u]nters [i]terator
};

