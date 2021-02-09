// Created by Petr Karnakov on 09.02.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <map>
#include <string>

namespace parse {

// Replaces '{key}' in `fmt` with `map[key]`
//
// Example:
// fmt: "a={key0} b={key1}"
// map: {"key0":"1", "key1":"s"}
// Returns:
// "a=1 b=s"
std::string SubstituteTemplate(
    const std::string& fmt, const std::map<std::string, std::string>& map);

// Returns values found in text `txt` by template `fmt`.
//
// Example:
// fmt: "a={key0} b={key1}"
// txt: "a=1 b=s"
// Returns:
// {"key0":"1", "key1":"s"}
std::map<std::string, std::string> ParseTemplate(
    const std::string& fmt, const std::string& txt);

} // namespace parse
