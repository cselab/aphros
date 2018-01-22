template <typename U>
unsigned bsr(U x)
{
  unsigned k;
#if __i386__ && USEASM
  __asm__("bsr %1, %0" : "=r"(k) : "r"(x));
#else
//  for (y = x, k = 0; y != 1; y >>= 1, k++);
  k = 0;
  do k++; while (x >>= 1);
  k--;
#endif
  return k;
}

