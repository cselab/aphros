#include <thread>
#include <iostream>
#include <fstream>
#include <mutex>
#include <queue>
#include <cassert>
#include <cstring>

template <class T>
class Queue {
 public:
  Queue() = default;
  Queue(const Queue&) = delete;
  Queue(Queue&&) = delete;
  void push(const T& a) {
    std::lock_guard<std::mutex> lock(m_);
    q_.push(a);
  }
  void pop() {
    std::lock_guard<std::mutex> lock(m_);
    return q_.pop();
  }
  bool empty() const { return q_.empty(); }
  T front() const { return q_.front(); }
  T back() const { return q_.back(); }

 private:
  std::queue<T> q_;
  mutable std::mutex m_;
};

void Listen() {
  Queue<std::string> q;
  std::string name = "chpipe";
  bool exit = false;
  while (!exit) {
    std::ifstream in(name);
    if (!in.good()) {
       std::cerr << "Error code: " << std::strerror(errno) << std::endl;
       std::terminate();
    }
    while (!exit && in) {
      std::string line;
      std::getline(in, line);
      if (!in) break;
      q.push(line);
      while (!q.empty()) {
        std::cerr << q.front() << std::endl;
        if (q.front() == "exit") {
          exit = true;
        }
        q.pop();
      }
    }
  }
}

int main() {
  std::thread t(Listen);
  //t.detach();
  t.join();
}
