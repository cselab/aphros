#ifndef RC_MODEL_H
#define RC_MODEL_H

class RCmodel {
public:
  RCmodel(unsigned symbols) : symbols(symbols) {}
  virtual ~RCmodel() {}

  // get frequency r for a symbol s and cumulative frequency l
  // of all symbols t < s
  virtual void encode(unsigned s, unsigned& l, unsigned& r) = 0;

  // get symbol frequency r and return symbol corresponding to
  // cumulative frequency l
  virtual unsigned decode(unsigned& l, unsigned& r) = 0;

  // divide range r by sum of all frequency counts
  virtual void normalize(unsigned &r) = 0;

  const unsigned symbols; // number of symbols
};

#endif
