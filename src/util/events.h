// Created by Petr Karnakov on 08.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <functional>
#include <map>
#include <string>

#include "parse/vars.h"

class Events {
 public:
  // Event handler. Takes <arg>.
  using Handler = std::function<void(std::string)>;

  // var: parameters to read events from, and apply events to
  // isroot: root block
  // islead: lead block
  Events(Vars& var, bool isroot, bool islead, bool verbose);
  // Parse events from var.String and put to events_
  void Parse();
  // Execute events due and remove from events_
  void Execute(double t);
  void AddHandler(std::string cmd, Handler);

 private:
  struct Event {
    double t;
    std::string cmd;
    std::string arg;
  };
  Vars& var_;
  const bool isroot_;
  const bool islead_;
  bool verbose_;
  std::map<std::string, Event> events_;
  // Mapping from command name <cmd> to handler.
  std::map<std::string, Handler> handlers_;
};
