// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <list>
#include <memory>
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
  class Sem { // [sem]aphore
   public:
    // Constructor
    // Advance list iterator, add new counter if needed, reset counter
    Sem(Suspender& owner_, std::string name = "");
    Sem(Sem&) = delete;
    Sem& operator=(Sem&) = delete;
    Sem(Sem&&) = default;
    // Destructor
    // If all nested stages done, next stage.
    // If all stages on current level done, remove current level
    ~Sem();
    // Next() without nested calls
    bool operator()(std::string suff = "" /*name suffix*/);
    // Next() with nested calls
    bool Nested(std::string suff = "" /*name suffix*/);
    std::string GetName() const {
      return name_;
    }
    void LoopBegin();
    void LoopBreak();
    void LoopEnd();
    template <class T>
    T* Get() {
      State& s = *owner_.pos_;
      std::unique_ptr<BaseHolder>& h = s.context;
      if (!h) {
        h = std::unique_ptr<BaseHolder>(new Holder<T>(new T()));
      }
      return dynamic_cast<Holder<T>*>(h.get())->Get();
    }
    template <class T>
    T* Get(T*) {
      return Get<T>();
    }
    template <class T>
    explicit operator T*() {
      return Get<T>();
    }

   private:
    Suspender& owner_;
    std::string name_;
    // Returns true if current stage needs execution
    // and advances stage counter
    bool Next(std::string suff = "" /*name suffix*/);
  };
  friend Sem;
  // Intializes list with auxiliary counter (-1,-1), sets iterator to it
  Suspender();
  Sem GetSem(std::string name = "");
  // Returns name+suff of current stage
  std::string GetCurName() const;
  // Converts counter list to string
  size_t GetDepth() const {
    return depth_;
  }
  std::string Print() const;
  // Returns true if there are unfinished levels
  bool Pending() const;

 private:
  class BaseHolder {
   public:
    virtual ~BaseHolder() {}
  };
  // Container for user-defined context object
  template <class T>
  class Holder : public BaseHolder {
   public:
    Holder(T* p) : p_(std::unique_ptr<T>(p)) {}
    T* Get() {
      return p_.get();
    }

   private:
    std::unique_ptr<T> p_;
  };
  struct State { // state of a suspended function
    int current; // current stage index
    int target; // target stage index
    int loop_begin = -1; // index of stage with LoopBegin
    int loop_end = -1; // index of stage with LoopEnd
    std::unique_ptr<BaseHolder> context;
    State(int current_, int target_) : current(current_), target(target_) {}
  };
  std::list<State> states_; // states of all suspended functions in nested call
  std::list<State>::iterator pos_; // current function
  bool nest_; // allow nested calls
  std::string curname_; // name+suff of current stage
  size_t depth_;
};
