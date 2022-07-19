// Created by Petr Karnakov on 09.02.2021
// Copyright 2021 ETH Zurich

#include "template.h"
#include "util/logger.h"

namespace parse {

std::string SubstituteTemplate(
    const std::string& fmt, const std::map<std::string, std::string>& map) {
  std::string res;
  size_t i = 0;
  std::string key;
  while (i < fmt.size()) {
    if (fmt[i] == '{') {
      key = "";
      ++i;
      for (; i < fmt.size(); ++i) {
        if (fmt[i] == '}') {
          break;
        }
        key += fmt[i];
      }
      ++i;
      fassert(map.count(key), "Key not found: " + key);
      res += map.at(key);
    } else {
      res.push_back(fmt[i]);
      ++i;
    }
  }
  return res;
}

std::map<std::string, std::string> ParseTemplate(
    const std::string& fmt, const std::string& txt) {
  std::map<std::string, std::string> res;
  size_t i = 0;
  size_t j = 0;
  std::string key;
  std::string value;
  auto is_ws = [](char c) { //
    return c == ' ' || c == '\n';
  };
  while (i < fmt.size() && j < txt.size()) {
    while (i < fmt.size() && is_ws(fmt[i])) {
      ++i;
    }
    while (j < txt.size() && is_ws(txt[j])) {
      ++j;
    }
    if (fmt[i] == '{') {
      if (fmt[i] == txt[j]) {
        std::stringstream msg;
        msg << "Expected different characters, got '" << fmt[i]
            << "' at position " << i << " of template and '" << txt[j]
            << "' at position " << j << " of input";
        fassert(false, msg.str());
      }
      key = "";
      value = "";
      ++i;
      for (; i < fmt.size(); ++i) {
        if (fmt[i] == '}') {
          break;
        }
        key += fmt[i];
      }
      ++i;
      while (i < fmt.size() && j < txt.size() && fmt[i] != txt[j]) {
        value += txt[j];
        ++j;
      }
      res[key] = value;
    } else {
      if (fmt[i] != txt[j]) {
        std::stringstream msg;
        msg << "Expected equal characters, got '" << fmt[i] << "' at position "
            << i << " of template and '" << txt[j] << "' at position " << j
            << " of input\n";
        fassert(false, msg.str());
      }
      ++i;
      ++j;
    }
  }
  return res;
}
} // namespace parse
