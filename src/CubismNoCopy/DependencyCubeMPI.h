/*
 *  DependencyCubeMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

struct Region {
  int s[3], e[3];

  bool operator<(const Region& r) const {
    const int* const a = &this->s[0];
    const int* const b = &r.s[0];

    for (int i = 0; i < 6; ++i)
      if (a[i] < b[i])
        return true;
      else if (a[i] > b[i])
        return false;

    return false;
  }

  void print() {
    printf("[%d-%d]x[%d-%d]x[%d-%d]", s[0], e[0], s[1], e[1], s[2], e[2]);
  };
};

template <typename Request>
class DependencyCubeMPI {
  enum State { PENDING = 0, READY = 1, USED = 2 };
  enum Side { MINUS = 0, PLUS = 1 };
  enum Direction { X = 0, Y = 1, Z = 2 };

  int n[3];
  State recvfaces[3][2], recvedges[3][2][2], recvcorners[2][2][2]; // events
  Request reqfaces[3][2], reqedges[3][2][2], reqcorners[2][2][2];

 public:
  set<Request> all_pending;

  struct Object {
    Region getRegion() {
      return Region();
    }
  };

  struct Face : Object {
    Direction d;
    Side s;

    Face() : d(X), s(MINUS) {}

    Face(Direction d, Side s) : d(d), s(s) {}

    Face& operator=(const Face& f) {
      d = f.d;
      s = f.s;
      return *this;
    }

    Region getRegion(const int n[3]) {
      Region r;

      const int d1 = (d + 1) % 3;
      const int d2 = (d + 2) % 3;

      r.s[d] = s * (n[d] - 1);
      r.s[d1] = 1;
      r.s[d2] = 1;

      r.e[d] = r.s[d] + 1;
      r.e[d1] = n[d1] - 1;
      r.e[d2] = n[d2] - 1;

      return r;
    }
  };

  struct Edge : Object {
    Direction d;
    Side a, b;

    Edge() : d(X), a(MINUS), b(MINUS) {}

    Edge(Direction d, Side a, Side b) : d(d), a(a), b(b) {}

    Edge& operator=(const Edge& f) {
      d = f.d;
      a = f.a;
      b = f.b;
      return *this;
    }

    Region getRegion(const int n[3]) {
      Region r;

      const int d1 = (d + 1) % 3;
      const int d2 = (d + 2) % 3;

      r.s[d] = 1;
      r.s[d1] = a * (n[d1] - 1);
      r.s[d2] = b * (n[d2] - 1);

      r.e[d] = n[d] - 1;
      r.e[d1] = r.s[d1] + 1;
      r.e[d2] = r.s[d2] + 1;

      return r;
    }
  };

  struct Corner : Object {
    Side x, y, z;

    Corner() : x(MINUS), y(MINUS), z(MINUS) {}
    Corner(Side x, Side y, Side z) : x(x), y(y), z(z) {}

    Corner& operator=(const Corner& f) {
      x = f.x;
      y = f.y;
      z = f.z;
      return *this;
    }

    Region getRegion(const int n[3]) {
      Region r;

      r.s[0] = x * (n[0] - 1);
      r.s[1] = y * (n[1] - 1);
      r.s[2] = z * (n[2] - 1);

      r.e[0] = r.s[0] + 1;
      r.e[1] = r.s[1] + 1;
      r.e[2] = r.s[2] + 1;

      return r;
    }
  };

  map<Region, set<Request> > pending;

  bool finalized;

 public:
  DependencyCubeMPI(const int nx, const int ny, const int nz) {
    n[0] = nx;
    n[1] = ny;
    n[2] = nz;
  }

  void face(Request req, int d, int s) {
    assert(finalized == false);

    reqfaces[d][s] = req;
    recvfaces[d][s] = PENDING;

    all_pending.insert(req);
  }

  void edge(Request req, int d, int a, int b) {
    assert(finalized == false);

    reqedges[d][b][a] = req;
    recvedges[d][b][a] = PENDING;

    all_pending.insert(req);
  }

  void corner(Request req, int x, int y, int z) {
    assert(finalized == false);

    reqcorners[z][y][x] = req;
    recvcorners[z][y][x] = PENDING;

    all_pending.insert(req);
  }

  void inspect() {
    for (typename map<Region, set<Request> >::iterator it = pending.begin();
         it != pending.end(); ++it) {
      int vol = 1;
      for (int i = 0; i < 3; ++i)
        vol *= it->first.e[i] - it->first.s[i];

      printf(
          "(%d) region %d %d %d, %d %d %d-->%d\n", vol, it->first.s[0],
          it->first.s[1], it->first.s[2], it->first.e[0], it->first.e[1],
          it->first.e[2], it->second.size());
    }
  }

  void make_dependencies(const bool /*isroot*/) {
    // set the center
    {
      Region r;

      r.s[0] = 1;
      r.s[1] = 1;
      r.s[2] = 1;

      r.e[0] = n[0] - 1;
      r.e[1] = n[1] - 1;
      r.e[2] = n[2] - 1;

      const bool nonemptyx = r.s[0] < r.e[0];
      const bool nonemptyy = r.s[1] < r.e[1];
      const bool nonemptyz = r.s[2] < r.e[2];

      if (nonemptyx && nonemptyy && nonemptyz) pending[r] = set<Request>();
    }

    for (int d = 0; d < 3; ++d)
      for (int s = 0; s < 2; ++s) {
        Face f((Direction)d, (Side)s);
        Region r = f.getRegion(n);

        pending[r] = set<Request>();

        if (recvfaces[d][s] == PENDING) pending[r].insert(reqfaces[d][s]);
      }

    for (int d = 0; d < 3; ++d)
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a) {
          Edge e((Direction)d, (Side)a, (Side)b);
          Region r = e.getRegion(n);

          pending[r] = set<Request>();

          if (recvedges[d][b][a] == PENDING)
            pending[r].insert(reqedges[d][b][a]);

          const int d1 = (d + 1) % 3;
          const int d2 = (d + 2) % 3;

          if (recvfaces[d1][a] == PENDING) pending[r].insert(reqfaces[d1][a]);

          if (recvfaces[d2][b] == PENDING) pending[r].insert(reqfaces[d2][b]);
        }

    for (int z = 0; z < 2; ++z)
      for (int y = 0; y < 2; ++y)
        for (int x = 0; x < 2; ++x) {
          Corner c((Side)x, (Side)y, (Side)z);
          Region r = c.getRegion(n);

          pending[r] = set<Request>();

          if (recvcorners[z][y][x] == PENDING)
            pending[r].insert(reqcorners[z][y][x]);

          if (recvfaces[0][x] == PENDING) pending[r].insert(reqfaces[0][x]);

          if (recvfaces[1][y] == PENDING) pending[r].insert(reqfaces[1][y]);

          if (recvfaces[2][z] == PENDING) pending[r].insert(reqfaces[2][z]);

          if (recvedges[0][z][y] == PENDING)
            pending[r].insert(reqedges[0][z][y]);

          if (recvedges[1][x][z] == PENDING)
            pending[r].insert(reqedges[1][x][z]);

          if (recvedges[2][y][x] == PENDING)
            pending[r].insert(reqedges[2][y][x]);
        }

    finalized = true;
  }

  void received(Request req) {
    assert(finalized == true);

    assert(all_pending.find(req) != all_pending.end());

    for (typename map<Region, set<Request> >::iterator it = pending.begin();
         it != pending.end(); ++it)
      if (it->second.find(req) !=
          it->second.end()) // this region needs this request
        it->second.erase(req);
  }

  vector<Region> avail() {
    assert(finalized == true);

    vector<Region> result;

    for (typename map<Region, set<Request> >::iterator it = pending.begin();
         it != pending.end(); ++it)
      if (it->second.size() == 0) // ready to be released
        result.push_back(it->first);

    for (vector<Region>::const_iterator it = result.begin(); it != result.end();
         ++it)
      pending.erase(*it);

    return result;
  }

  int pendingcount() const {
    return pending.size();
  }

  void prepare() {
    finalized = false;

    for (int d = 0; d < 3; ++d)
      for (int s = 0; s < 2; ++s)
        recvfaces[d][s] = READY;

    for (int d = 0; d < 3; ++d)
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          recvedges[d][b][a] = READY;

    for (int z = 0; z < 2; ++z)
      for (int y = 0; y < 2; ++y)
        for (int x = 0; x < 2; ++x)
          recvcorners[z][y][x] = READY;
  }
};
