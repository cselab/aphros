#ifndef RC_QSMODEL_H
#define RC_QSMODEL_H

#include "rcmodel.h"

class RCqsmodel : public RCmodel {
public:
  // initialization of model
  // compress: true for compression, false for decompression
  // symbols:  number of symbols
  // bits:     log2 of total frequency count (must be <= 16)
  // period:   max symbols between normalizations (must be < 1<<(bits+1))
  RCqsmodel(bool compress, unsigned symbols, unsigned bits = 16, unsigned period = 0x400);
  ~RCqsmodel();

  // reinitialize model
  void reset();

  // get frequencies for a symbol s
  void encode(unsigned s, unsigned& l, unsigned& r);

  // return symbol corresponding to cumulative frequency l
  unsigned decode(unsigned& l, unsigned& r);

  //
  void normalize(unsigned &r);

private:

  void update();
  void update(unsigned s);

  const unsigned  bits;    // number of bits of precision for frequencies

  unsigned  left;          // number of symbols until next normalization
  unsigned  more;          // number of symbols with larger increment
  unsigned  incr;          // increment per update
  unsigned  rescale;       // current interval between normalizations
  unsigned  targetrescale; // target interval between rescales
  unsigned* symf;          // array of partially updated frequencies
  unsigned* cumf;          // array of cumulative frequencies

  unsigned  searchshift;   // difference of frequency bits and table bits
  unsigned* search;        // structure for searching on decompression
};

#include "rcqsmodel.inl"

#endif
