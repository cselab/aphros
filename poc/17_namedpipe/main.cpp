// Created by Petr Karnakov on 25.10.2019
// Copyright 2019 ETH Zurich

#include <atomic>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>

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
  bool empty() const {
    return q_.empty();
  }
  T front() const {
    return q_.front();
  }
  T back() const {
    return q_.back();
  }

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
    std::this_thread::yield();
    // std::cerr << "opening pipe '" << name << "'" << std::endl;
    std::ifstream in(name);
    if (!in.good()) {
      Abort();
    }
    while (!flagexit && in) {
      std::string line;
      std::getline(in, line);
      if (!in) break;
      queue.push(line);
      std::this_thread::yield();
    }
  }
}

void Send() {
  static double base = 0;
  std::string name = "out";
  while (!flagexit) {
    while (queue.empty()) {
      std::this_thread::yield();
    }
    if (queue.front() == "exit") {
      flagexit = true;
    }

    std::ofstream out(name);
    if (!out.good()) {
      Abort();
    }

    std::cerr << "send: " << queue.front() << std::endl;
    std::stringstream s(queue.front());
    std::string cmd;
    s >> cmd;
    if (cmd == "sum") {
      out << cmd << std::endl;
      double a, b;
      s >> a >> b;
      out << a + b << std::endl;
    } else if (cmd == "field") {
      out << cmd << std::endl;
      int nx = 10;
      int ny = 10;
      int nz = 1;
      out << nx << " " << ny << " " << nz << std::endl;
      for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
          double a = std::sin(x * 0.5 + base) * std::sin(y * 0.5 + base * 0.5);
          out << a << " ";
        }
      }
      out << std::endl;
      out << std::endl;
    } else if (cmd == "step") {
      base += 0.1;
      out << cmd << std::endl;
      out << "base=" << base << std::endl;
    } else {
      out << queue.front() << std::endl;
    }
    queue.pop();
  }
}

int main() {
  std::thread t1(Listen);
  std::thread t2(Send);
  // t.detach();
  t1.join();
  t2.join();
}
