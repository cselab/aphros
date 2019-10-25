#include <thread>
#include <iostream>
#include <fstream>
#include <mutex>
#include <queue>
#include <cassert>
#include <cstring>
#include <atomic>

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

Queue<std::string> queue;

std::atomic<bool> flagexit(false);

void Abort() {
  std::cerr << "Error code: " << std::strerror(errno) << std::endl;
  std::terminate();
}

void Listen() {
  std::string name = "in";
  while (!flagexit) {
    std::cerr << "opening pipe '" << name << "'" << std::endl;
    std::ifstream in(name);
    if (!in.good()) { Abort(); }
    while (!flagexit && in) {
      std::string line;
      std::getline(in, line);
      if (!in) break;
      queue.push(line);
    }
  }
}

void Send() {
  std::string name = "out";
  while (!flagexit) {
    std::cerr << "opening pipe '" << name << "'" << std::endl;
    std::ofstream out(name);
    if (!out.good()) { Abort(); }
    while (out && !flagexit) {
      while (queue.empty()) {
        std::this_thread::yield();
      }
      std::cerr << "send " << queue.front() << std::endl;
      if (queue.front() == "exit") {
        flagexit = true;
      }
      out << queue.front() << std::endl;
      queue.pop();
    }
  }
}

int main() {
  std::thread t1(Listen);
  std::thread t2(Send);
  //t.detach();
  t1.join();
  t2.join();
}
