/*
 *  StencilInfo.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 11/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cassert>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

using namespace std;

struct StencilInfo {
  int sx, sy, sz, ex, ey, ez;
  vector<int> selcomponents;

  bool tensorial;

  StencilInfo() {}

  StencilInfo(
      int sx, int sy, int sz, int ex, int ey, int ez, bool tensorial,
      const int ncomponents, ...)
      : sx(sx)
      , sy(sy)
      , sz(sz)
      , ex(ex)
      , ey(ey)
      , ez(ez)
      , selcomponents()
      , tensorial(tensorial) {
    assert(ncomponents > 0);

    va_list ev_list;

    va_start(ev_list, ncomponents);
    for (int i = 0; i < ncomponents; i++)
      selcomponents.push_back(va_arg(ev_list, int));
    va_end(ev_list);

    if (!isvalid()) {
      cout << "Stencilinfo instance not valid. Aborting\n";
      abort();
    }
  }

  StencilInfo(const StencilInfo& c)
      : sx(c.sx)
      , sy(c.sy)
      , sz(c.sz)
      , ex(c.ex)
      , ey(c.ey)
      , ez(c.ez)
      , selcomponents(c.selcomponents)
      , tensorial(c.tensorial) {}

  // serialize all fields as vector<int>
  vector<int> _all() const {
    int extra[] = {sx, sy, sz, ex, ey, ez, (int)tensorial};
    vector<int> all(selcomponents);
    all.insert(all.end(), extra, extra + sizeof(extra) / sizeof(int));

    return all;
  }

  // lexical order by serialized fields _all()
  bool operator<(StencilInfo s) const {
    vector<int> me = _all(), you = s._all();

    const int N = min(me.size(), you.size());

    for (int i = 0; i < N; ++i)
      if (me[i] < you[i])
        return true;
      else if (me[i] > you[i])
        return false;

    return me.size() < you.size();
  }

  bool isvalid() const {
    const bool not0 = selcomponents.size() == 0;
    const bool not1 = sx > 0 || ex <= 0 || sx > ex;
    const bool not2 = sy > 0 || ey <= 0 || sy > ey;
    const bool not3 = sz > 0 || ez <= 0 || sz > ez;

    return !(not0 || not1 || not2 || not3);
  }
};
