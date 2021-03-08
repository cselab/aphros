// Created by Petr Karnakov on 08.03.2021
// Copyright 2021 ETH Zurich

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <ios>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "util/logger.h"

struct CodeBlock {
  std::string name;
  std::string content;
};

static std::string Strip(std::string s) {
  size_t b = 0;
  size_t e = s.size();
  while (b < e && std::isspace(s[b])) {
    ++b;
  }
  while (b < e && std::isspace(s[e - 1])) {
    --e;
  }
  return s.substr(b, e - b);
}

// Extacts the next code block from stream f.
// Returns true if found.
auto ParseCodeBlock(std::istream& fin) {
  struct Res {
    bool found;
    CodeBlock block;
    std::string error;
  };
  const auto flags = fin.flags();
  fin >> std::noskipws;
  CodeBlock block;
  enum class S { begin, name, content, exit };
  S s = S::begin; // state
  char c = ' ';
  int braces = 0;
  std::string error;
  while (fin && s != S::exit && error.empty()) {
    auto next = [&fin, &c]() { fin >> c; };
    switch (s) {
      case S::begin: {
        if (std::isspace(c)) {
          next();
        } else {
          s = S::name;
        }
        break;
      }
      case S::name: {
        if (c == '{') {
          s = S::content;
          next();
        } else if (c == '}') {
          error = "name cannot contain '}', unmatched brace?";
        } else {
          block.name += c;
          next();
        }
        break;
      }
      case S::content: {
        if (c == '}' && braces == 0) {
          s = S::exit;
          next();
        } else {
          block.content += c;
          if (c == '{') {
            ++braces;
          } else if (c == '}') {
            --braces;
          }
          next();
        }
        break;
      }
      case S::exit: {
        break;
      }
    }
  }
  if (error.empty() && s != S::exit &&
      (block.name.size() || block.content.size())) {
    error = "unexpected end of stream";
  }
  if (!error.empty()) {
    error = std::string(__func__) + ": " + error +
            "\nstopping at block:\nname='" + Strip(block.name) + "' content='" +
            Strip(block.content) + "'";
  }
  fin.flags(flags);
  block.name = Strip(block.name);
  block.content = Strip(block.content);
  return Res{s == S::exit, block, error};
}

// Parses a stream and returns a list of code blocks.
std::vector<CodeBlock> ParseCodeBlocks(std::istream& f) {
  std::vector<CodeBlock> blocks;
  std::string error;
  while (true) {
    auto res = ParseCodeBlock(f);
    if (!res.error.empty()) {
      error = res.error;
      break;
    }
    if (!res.found) {
      break;
    }
    blocks.push_back(res.block);
  }

  if (!error.empty()) {
    std::string o;
    o += __func__;
    o += ": error after parsing " + std::to_string(blocks.size()) + " blocks\n";
    o += error;
    if (blocks.size()) {
      o += "\nlast successful block:\n";
      o += "name='" + Strip(blocks.back().name) + "' ";
      o += "content='" + Strip(blocks.back().content) + "'";
    }
    fassert(false, o);
  }
  return blocks;
}
