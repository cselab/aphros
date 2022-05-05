// Created by Petr Karnakov on 05.05.2022
// Copyright 2022 ETH Zurich

#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

template <class T, class Iterator>
T EvalExpr(Iterator begin, Iterator end);

template <class T>
T EvalExpr(const std::string& expr);
