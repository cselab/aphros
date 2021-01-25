// Created by Petr Karnakov on 08.04.2019
// Copyright 2019 ETH Zurich

#include <iomanip>
#include <iostream>
#include <sstream>

#include "events.h"
#include "parse/parser.h"

Events::Events(Vars& var, bool isroot, bool islead)
    : var_(var), isroot_(isroot), islead_(islead) {
  AddHandler("echo", [isroot = isroot_](std::string arg) {
    if (isroot) {
      std::cout << arg << std::endl;
    }
  });
  AddHandler("set", [this](std::string arg) {
    Parser p(var_);
    if (islead_) {
      p.Run("set " + arg);
    }
  });
}

void Events::AddHandler(std::string cmd, Handler h) {
  handlers_.emplace(cmd, h);
}

// Parses events in var.String selecting variables by names:
// evN <time> <cmd> <arg>
// Example:
// ev0 0.1 set double prelax 1
void Events::Parse() {
  // Check at least first nmax indices and all contiguous
  int n = 0;
  const int nmax = 100;
  while (true) {
    std::string k = "ev" + std::to_string(n);
    if (auto p = var_.String.Find(k)) {
      Event e;
      std::stringstream b(*p); // buf
      b >> std::skipws;

      // format: <time> <cmd> <arg>
      b >> e.t;
      b >> e.cmd;

      char c;
      // Read first non-ws character
      b >> c;
      // Read remaining line
      std::string s;
      std::getline(b, s);
      e.arg = c + s;

      events_.emplace(k, e);
    } else if (n > nmax) {
      break;
    }
    ++n;
  }

  if (isroot_) {
    std::cout << "Found events: \n=====" << std::endl;
    for (auto p : events_) {
      Event& e = p.second;
      std::cout << p.first << " " << e.t << " " << e.cmd << " " << e.arg
                << std::endl;
    }
    std::cout << "=====" << std::endl;
  }
}

// Executes events in events_
// set <type> <key> <value>
// echo <string>
void Events::Execute(double t) {
  for (auto it = events_.begin(); it != events_.end();) {
    auto& event = it->second;
    const std::string cmd = event.cmd;
    const std::string arg = event.arg;

    if (t >= event.t) {
      if (isroot_) {
        std::cout << std::fixed << std::setprecision(8)
                  << "Event at t=" << event.t << ": " << cmd << " " << arg
                  << std::endl;
      }
      auto ith = handlers_.find(cmd);
      if (ith != handlers_.end()) {
        ith->second(arg);
      } else {
        throw std::runtime_error("ExecEvents(): Unknown command '" + cmd + "'");
      }
      it = events_.erase(it);
    } else {
      ++it;
    }
  }
}
