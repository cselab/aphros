#pragma once

#include <string>
#include <utility>


class Timer {
 public:
  Timer(std::string name, double timeout /*sec*/, size_t batch);
  Timer(std::string name, double timeout);
  Timer(std::string name);
  ~Timer() {}
  std::string GetName() const;
  // Runs F() until reaching the timeout in batches of b_ iterations.
  // Returns execution time per iteration and number of iterations.
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
