// Created by Petr Karnakov on 20.02.2020
// Copyright 2020 ETH Zurich

#pragma once

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
std::pair<bool, CodeBlock> ParseCodeBlock(std::istream& fin) {
  const auto flags = fin.flags();
  fin >> std::noskipws;
  CodeBlock block;
  enum class S { begin, name, content, exit };
  S s = S::begin; // state
  char c = ' ';
  int braces = 0;
  try {
    while (fin && s != S::exit) {
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
            throw std::runtime_error(
                "name cannot contain '}', unmatched brace?");
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
    if (s != S::exit && (block.name.size() || block.content.size())) {
      throw std::runtime_error("unexpected end of stream");
    }
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(
        std::string(__func__) + ": " + e.what() +
        "\nstopping at block:\nname='" + Strip(block.name) + "' content='" +
        Strip(block.content) + "'");
  }
  fin.flags(flags);
  block.name = Strip(block.name);
  block.content = Strip(block.content);
  return {s == S::exit, block};
}

// Parses a stream and returns a list of code blocks.
std::vector<CodeBlock> ParseCodeBlocks(std::istream& f) {
  std::vector<CodeBlock> blocks;
  try {
    while (true) {
      auto pair = ParseCodeBlock(f);
      if (!pair.first) {
        break;
      }
      blocks.push_back(pair.second);
    }
  } catch (const std::runtime_error& e) {
    std::string o;
    o += __func__;
    o += ": error after parsing " + std::to_string(blocks.size()) + " blocks\n";
    o += e.what();
    if (blocks.size()) {
      o += "\nlast successful block:\n";
      o += "name='" + Strip(blocks.back().name) + "' ";
      o += "content='" + Strip(blocks.back().content) + "'";
    }
    throw std::runtime_error(o);
  }
  return blocks;
}
