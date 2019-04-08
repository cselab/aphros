#pragma once

#include <list>
#include <string>

// TODO: Sequence like 
//   GetSem();
//   sem = GetSem();
//   sem();
// causes an infinite loop
// (see todo in tests/suspender)
//
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
    int lb; // loop begin
    int le; // loop end
    U(int c, int t) : c(c), t(t), lb(-1), le(-1) {}
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
    // Next() without nested calls
    bool operator()(std::string suff="" /*name suffix*/);
    // Next() with nested calls
    bool Nested(std::string suff="" /*name suffix*/);
    std::string GetName() const { return name_; }
    void LoopBegin();
    void LoopBreak();
    void LoopEnd();
   private:
    Suspender& p; // parent
    std::string name_;
    // Returns true if current stage needs execution
    // and advances stage counter
    bool Next(std::string suff="" /*name suffix*/);
  };
  friend Sem;
  // Intializes list with auxiliary counter (-1,-1), sets iterator to it
  Suspender();
  Sem GetSem(std::string name="");
  // Returns name+suff of current stage
  std::string GetCurName() const;
  // Converts counter list to string
  size_t GetDepth() const { return depth_; }
  std::string Print() const;
  // Returns true if there are unfinished levels 
  bool Pending() const;

 private:
  using LU = std::list<U>;
  LU lu_;      // [l]ist of co[u]nters
  LU::iterator lui_; // [l]ist of co[u]nters [i]terator
  bool nest_; // allow nested calls
  std::string curname_; // name+suff of current stage
  size_t depth_;
};
