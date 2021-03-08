// Created by Petr Karnakov on 20.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

struct CodeBlock {
  std::string name;
  std::string content;
};

// Parses a stream and returns a list of code blocks.
std::vector<CodeBlock> ParseCodeBlocks(std::istream& f);
