#include <cassert>
#include "rcqsmodel.h"

// table size for binary search
#define TBLSHIFT 7

RCqsmodel::RCqsmodel(bool compress, unsigned symbols, unsigned bits, unsigned period) : RCmodel(symbols), bits(bits), targetrescale(period)
{
  assert(bits <= 16);
  assert(targetrescale < (1u << (bits + 1)));

  unsigned n = symbols;
  symf = new unsigned[n + 1];
  cumf = new unsigned[n + 1];
  cumf[0] = 0;
  cumf[n] = 1 << bits;
  if (compress)
    search = 0;
  else {
    searchshift = bits - TBLSHIFT;
    search = new unsigned[(1 << TBLSHIFT) + 1];
  }
  reset();
}

RCqsmodel::~RCqsmodel()
{
  delete [] symf;
  delete [] cumf;
  delete [] search;
}

// reinitialize model
void RCqsmodel::reset()
{
  unsigned n = symbols;
  rescale = (n >> 4) | 2;
  more = 0;
  unsigned f = cumf[n] / n;
  unsigned m = cumf[n] % n;
  for (unsigned i = 0; i < m; i++)
    symf[i] = f + 1;
  for (unsigned i = m; i < n; i++)
    symf[i] = f;
  update();
}

// return symbol corresponding to cumulative frequency l
unsigned RCqsmodel::decode(unsigned& l, unsigned& r)
{
  unsigned i = l >> searchshift;
  unsigned s = search[i];
  unsigned h = search[i + 1] + 1;

  // find symbol via binary search
  while (s + 1 < h) {
    unsigned m = (s + h) / 2;
    if (l < cumf[m])
      h = m;
    else
      s = m;
  }

  l = cumf[s];
  r = cumf[s + 1] - l;
  update(s);

  return s;
}

// update probability table
void RCqsmodel::update()
{
  if (more) { // we have some more symbols before update
    left = more;
    more = 0;
    incr++;
    return;
  }
  if (rescale != targetrescale) {
    rescale *= 2;
    if (rescale > targetrescale)
      rescale = targetrescale;
  }

  // update symbol frequencies
  unsigned n = symbols;
  unsigned cf = cumf[n];
  unsigned count = cf;
  for (unsigned i = n; i--; ) {
    unsigned sf = symf[i];
    cf -= sf;
    cumf[i] = cf;
    sf = (sf >> 1) | 1;
    count -= sf;
    symf[i] = sf;
  }
  // assert(cf == 0);

  // count is now difference between target cumf[n] and sum of symf;
  // next actual rescale happens when sum of symf equals cumf[n]
  incr = count / rescale;
  more = count % rescale;
  left = rescale - more;

  // build lookup table for fast symbol searches
  if (search)
    for (unsigned i = n, h = 1 << TBLSHIFT; i--; h = cumf[i] >> searchshift)
      for (unsigned l = cumf[i] >> searchshift; l <= h; l++)
        search[l] = i;
}

// update frequency for symbol s
void RCqsmodel::update(unsigned s)
{
  if (!left)
    update();
  left--;
  symf[s] += incr;
}
