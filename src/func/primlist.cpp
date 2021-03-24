// Created by Petr Karnakov on 02.12.2020
// Copyright 2020 ETH Zurich

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "util/format.h"
#include "util/logger.h"

#include "primlist.ipp"

#define X(dim) template struct UPrimList<generic::Vect<double, dim>>;
MULTIDIMX
