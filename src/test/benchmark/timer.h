#pragma once

#include <string>
#include <utility>

class Timer {
 public:
  // name: returned by GetName()
  // timeout: time in seconds for which to repeat F() 
  // batch: number of calls F() at every iteration
  Timer(std::string name, double timeout, size_t batch);
  // batch=1
  Timer(std::string name, double timeout);
  // timeout=0.01, batch=1
  Timer(std::string name);
  virtual ~Timer() {}
  std::string GetName() const;
  // Runs F() until reaching the timeout in batches of b_ calls.
  // Returns execution time per call and number of calls.
  std::pair<double, size_t> Run();

 private:
  // Function to evaluate.
  // Implementation note: use a volatile variable to prevent optimization.
  virtual void F() = 0;
  void B();

  std::string n_; // name
  double to_;     // timeout
  size_t b_;      // batch size
};
