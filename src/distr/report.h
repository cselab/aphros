// Created by Petr Karnakov on 21.07.2018
// Copyright 2018 ETH Zurich

#include <iosfwd>
#include <map>
#include <string>

// timings: map from stage name to timing, stage name composed from
//   arguments of Sem() separated by ' --> '
// Example:
//   "fluid --> step --> init" : 0.1
void ParseReport(
    const std::map<std::string, double>& timings, std::ostream& out);
