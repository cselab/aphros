// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <exception>
#include <string>

#include "util/logger.h"

class UnsteadySolver {
 public:
  UnsteadySolver(double t, double dt) : t_(t), dt_(dt) {}
  virtual ~UnsteadySolver() {}
  virtual void StartStep() {}
  virtual void FinishStep() {
    IncTime();
  }
  virtual double GetTime() const {
    return t_;
  }
  virtual void SetTime(double t) {
    t_ = t;
  }
  virtual double GetTimeStep() const {
    return dt_;
  }
  virtual void SetTimeStep(double dt) {
    dt_ = dt;
  }

 protected:
  virtual void IncTime() {
    t_ += dt_;
  }

 private:
  double t_;
  double dt_;
};

class UnsteadyIterativeSolver : public UnsteadySolver {
 public:
  UnsteadyIterativeSolver(double t, double dt) : UnsteadySolver(t, dt), i_(0) {}
  virtual void MakeIteration() = 0;
  virtual double GetError() const {
    return 0.;
  }
  virtual size_t GetIter() const {
    return i_;
  }
  virtual void StartStep() override {
    ClearIter();
  }

 protected:
  void IncIter() {
    ++i_;
  }
  void ClearIter() {
    i_ = 0;
  }

 private:
  size_t i_;
};

// Between StartStep() and FinishStep():
//   iter_curr -- last iteration of next step
//   iter_prev -- previous iteration of next step
//   time_curr -- current step
//   time_prev -- previous step
// After FinishStep() call:
//   iter_curr, iter_prev -- undefined
//   time_curr -- current step
//   time_prev -- previous step
enum class Step { time_curr, time_prev, iter_curr, iter_prev };

template <class T>
struct StepData {
  T time_curr, time_prev, iter_curr, iter_prev;
  const T& Get(Step l) const {
    switch (l) {
      case Step::time_curr: {
        return time_curr;
      }
      case Step::time_prev: {
        return time_prev;
      }
      case Step::iter_curr: {
        return iter_curr;
      }
      case Step::iter_prev: {
        return iter_prev;
      }
      default: {
        fassert(false, "StepData::Get(): Unknown layer");
      }
    }
  }
  T& Get(Step l) {
    return const_cast<T&>(const_cast<const StepData*>(this)->Get(l));
  }
};

std::string GetName(Step);

// Convection scheme
enum class ConvSc { fou, cd, sou, quick, superbee };
// Convection scheme name
std::string GetName(ConvSc sc);
// Convection scheme by name.
ConvSc GetConvSc(std::string s);
