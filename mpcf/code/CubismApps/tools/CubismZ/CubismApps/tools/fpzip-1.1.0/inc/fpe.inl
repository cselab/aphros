template <typename T, typename U, unsigned mbits, unsigned ebits>
FPE<T, U, mbits, ebits>
FPE<T, U, mbits, ebits>::operator +(FPE<T, U, mbits, ebits> x) const
{
  // addition
  if (equalsign(x))
    return u > x.u ? add(x.u) : x.add(u);
  else {
    x = -x;
    return u > x.u ? sub(x.u) : u < x.u ? -x.sub(u) : FPE(U(0));
  }
}

template <typename T, typename U, unsigned mbits, unsigned ebits>
FPE<T, U, mbits, ebits>
FPE<T, U, mbits, ebits>::operator -(FPE<T, U, mbits, ebits> x) const
{
  // subtraction
  if (equalsign(x))
    return u > x.u ? sub(x.u) : u < x.u ? -x.sub(u) : FPE(U(0));
  else {
    x = -x;
    return u > x.u ? add(x.u) : x.add(u);
  }
}

template <typename T, typename U, unsigned mbits, unsigned ebits>
FPE<T, U, mbits, ebits>
FPE<T, U, mbits, ebits>::add(U g) const
{
  // effective addition f + g, |f| >= |g|
  U f = u;
  unsigned k = (f >> m) - (g >> m);
  if (k <= m) {
    U c = (e - (f & (e - 1)));
    U d = (e + (g & (e - 1))) >> k;
    if (c < d)
      d = (c + d) >> 1;
    f += d;
  }
  return FPE(f);
}

template <typename T, typename U, unsigned mbits, unsigned ebits>
FPE<T, U, mbits, ebits>
FPE<T, U, mbits, ebits>::sub(U g) const
{
  // effective subtraction f - g, |f| >= |g|
  U f = u;
  unsigned k = (f >> m) - (g >> m);
  if (k <= m) {
    U c = ((0 + (f & (e - 1))) << 1);
    U d = ((e + (g & (e - 1))) << 1) >> k;
    while (c < d) {
      d += d - c;
      c += e << 1;
    }
    f -= d >> 1;
  }
  return FPE(f);
}
