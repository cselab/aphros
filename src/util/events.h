#pragma once

#include <string>
#include <map>

#include "parse/vars.h"

class Events {
 public:
  // var: parameters to read events from, and apply events to
  // isroot: root block
  // islead: lead block
  Events(Vars& var, bool isroot, bool islead);
  // Parse events from var.String and put to ev_
  void Parse();
  // Exec events due and remove from ev_
  void Exec(double t);

 private:
  struct Event {
    double t;
    std::string cmd;
    std::string arg;
  };
  Vars& var_;
  const bool isroot_, islead_;
  std::map<std::string, Event> ev_;
};

