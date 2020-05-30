/*
 *  PUPKernelsMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#ifndef PUPKERNELSMPI_H_XBFRP5V0
#define PUPKERNELSMPI_H_XBFRP5V0

#include <cstring>

namespace PUPkernelsMPI {
inline void pack(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    int* selected_components, const int ncomponents, const int xstart,
    const int ystart, const int zstart, const int xend, const int yend,
    const int zend, const int bx, const int by) {
  for (int idst = 0, iz = zstart; iz < zend; ++iz)
    for (int iy = ystart; iy < yend; ++iy)
      for (int ix = xstart; ix < xend; ++ix) {
        const Real* src = srcbase + gptfloats * (ix + bx * (iy + by * iz));

        // bgq: s_c[ic] = ic! -> memcpy or stripes...
        for (int ic = 0; ic < ncomponents; ic++, idst++)
          dst[idst] = src[selected_components[ic]];
      }
}

inline void pack_soa(
    const Real* const __restrict srcbase, Real* const __restrict dstbase,
    const int xstart, const int ystart, const int zstart, const int xend,
    const int yend, const int zend, const int bx, const int by) {
  Real* dst = dstbase;
  for (int iz = zstart; iz < zend; ++iz) {
    for (int iy = ystart; iy < yend; ++iy) {
      const Real* src = srcbase + xstart + bx * (iy + by * iz);
      for (int ix = xstart; ix < xend; ++ix) {
        *dst++ = *src++;
      }
    }
  }
}

inline void pack_soa_ncomp(
    const Real* const __restrict srcbase, Real* const __restrict dstbase,
    const int ncomp, const int xstart, const int ystart, const int zstart,
    const int xend, const int yend, const int zend, const int bx,
    const int by) {
  Real* dst = dstbase;
  for (int iz = zstart; iz < zend; ++iz) {
    for (int iy = ystart; iy < yend; ++iy) {
      for (int ix = xstart; ix < xend; ++ix) {
        const Real* src = srcbase + ncomp * (ix + bx * (iy + by * iz));
        for (int c = 0; c < ncomp; ++c) {
          *dst++ = *src++;
        }
      }
    }
  }
}

inline void pack_stripes1(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    const int selstart, const int selend, const int xstart, const int ystart,
    const int zstart, const int xend, const int yend, const int zend,
    const int bx, const int by) {
  for (int idst = 0, iz = zstart; iz < zend; ++iz)
    for (int iy = ystart; iy < yend; ++iy)
      for (int ix = xstart; ix < xend; ++ix) {
        const Real* src = srcbase + gptfloats * (ix + bx * (iy + by * iz));

        for (int ic = selstart; ic < selend; ic++, idst++)
          dst[idst] = src[ic];
      }
}

inline void pack_stripes_(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    const int selstart, const int selend, const int xstart, const int ystart,
    const int zstart, const int xend, const int yend, const int zend,
    const int bx, const int by) {
  const int seldiff = selend - selstart;
  const int nbytes = seldiff * sizeof(Real);
  for (int idst = 0, iz = zstart; iz < zend; ++iz) {
    for (int iy = ystart; iy < yend; ++iy) {
      for (int ix = xstart; ix < xend; ++ix) {
        const Real* src = srcbase + gptfloats * (ix + bx * (iy + by * iz));

        std::memcpy(&dst[idst], &src[selstart], nbytes);
        idst += seldiff;
      }
    }
  }
}

inline void pack_stripes_x(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    const int selstart, const int selend, const int xstart, const int ystart,
    const int zstart, const int xend, const int yend, const int zend,
    const int bx, const int by) {
  const int seldiff = selend - selstart;
  const int nbytes = seldiff * sizeof(Real);
  const int _BS_XY_ = bx * by;
  int idst = 0;
  for (int iz = zstart; iz < zend; ++iz) {
    const int iz_off = _BS_XY_ * iz;
    for (int iy = ystart; iy < yend; ++iy) {
      const int iy_off = bx * iy;
      for (int ix = xstart; ix < xend; ++ix) {
        const Real* src = srcbase + gptfloats * (ix + iy_off + iz_off);
        std::memcpy(&dst[idst], &src[selstart], nbytes);
        idst += seldiff;
      }
    }
  }
}

#ifdef memcpy2
#undef memcpy2
#endif /* memcpy2 */

#include "QPXEMU.h"
#ifdef __bgq__
#include <builtins.h>
#define memcpy2(a, b, c) __bcopy((b), (a), (c))
#else
#define memcpy2(a, b, c) std::memcpy((a), (b), (c))
#endif

inline void pack_stripes_unroll0(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    const int selstart, const int selend, const int xstart, const int ystart,
    const int zstart, const int xend, const int yend, const int zend,
    const int bx, const int by) {
  const int seldiff = selend - selstart;
  const int nbytes = seldiff * sizeof(Real);
  for (int idst = 0, iz = zstart; iz < zend; ++iz) {
    for (int iy = ystart; iy < yend; ++iy) {
      int xentries = xend - xstart;
      // int unroll = 4;

      int repeat = (xentries / 4);
      int left = (xentries % 4);

      int ix = xstart;
      while (repeat-- > 0) {
        const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));
        const Real* src1 = src0 + gptfloats;
        const Real* src2 = src1 + gptfloats;
        const Real* src3 = src2 + gptfloats;

        memcpy2(
            (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
        memcpy2(
            (char*)&dst[idst + 1 * seldiff], (char*)&src1[selstart], nbytes);
        memcpy2(
            (char*)&dst[idst + 2 * seldiff], (char*)&src2[selstart], nbytes);
        memcpy2(
            (char*)&dst[idst + 3 * seldiff], (char*)&src3[selstart], nbytes);
        idst += 4 * seldiff;
        ix += 4;
      }

      if (left == 3) {
        const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));
        const Real* src1 = src0 + gptfloats;
        const Real* src2 = src1 + gptfloats;

        memcpy2(
            (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
        memcpy2(
            (char*)&dst[idst + 1 * seldiff], (char*)&src1[selstart], nbytes);
        memcpy2(
            (char*)&dst[idst + 2 * seldiff], (char*)&src2[selstart], nbytes);
        idst += 3 * seldiff;
        ix += 3;
      } else if (left == 2) {
        const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));
        const Real* src1 = src0 + gptfloats;

        memcpy2(
            (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
        memcpy2(
            (char*)&dst[idst + 1 * seldiff], (char*)&src1[selstart], nbytes);
        idst += 2 * seldiff;
        ix += 2;
      } else /* left == 1 */
      {
        const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));

        memcpy2(
            (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
        idst += 1 * seldiff;
        ix += 1;
      }
    }
  }
}

inline void pack_stripesxx(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    const int selstart, const int selend, const int xstart, const int ystart,
    const int zstart, const int xend, const int yend, const int zend,
    const int bx, const int by) {
  //	printf("xstart/end = (%d, %d)\n", xstart, xend);
  const int seldiff = selend - selstart;
  const int nbytes = seldiff * sizeof(Real);
  for (int idst = 0, iz = zstart; iz < zend; ++iz) {
    for (int iy = ystart; iy < yend; ++iy) {
      if ((xend - xstart) % 4 != 0) {
        for (int ix = xstart; ix < xend; ++ix) {
          const Real* src = srcbase + gptfloats * (ix + bx * (iy + by * iz));

          memcpy2((char*)&dst[idst], (char*)&src[selstart], nbytes);
          idst += seldiff;
        }
      }
      /*			else if ((xend - xstart)%8 == 0)
                              {
                                      for(int ix=xstart; ix<xend; ix+=8)
                                      {
                                              const Real * src0 = srcbase +
         gptfloats*(ix + bx*(iy + by*iz)); const Real * src1 = src0 + gptfloats;
                                              const Real * src2 = src1 +
         gptfloats; const Real * src3 = src2 + gptfloats; const Real * src4 =
         src3 + gptfloats; const Real * src5 = src4 + gptfloats; const Real *
         src6 = src5 + gptfloats; const Real * src7 = src6 + gptfloats;

                                              memcpy2((char
         *)&dst[idst+0*seldiff], (char *)&src0[selstart], nbytes); memcpy2((char
         *)&dst[idst+1*seldiff], (char *)&src1[selstart], nbytes); memcpy2((char
         *)&dst[idst+2*seldiff], (char *)&src2[selstart], nbytes); memcpy2((char
         *)&dst[idst+3*seldiff], (char *)&src3[selstart], nbytes); memcpy2((char
         *)&dst[idst+4*seldiff], (char *)&src4[selstart], nbytes); memcpy2((char
         *)&dst[idst+5*seldiff], (char *)&src5[selstart], nbytes); memcpy2((char
         *)&dst[idst+6*seldiff], (char *)&src6[selstart], nbytes); memcpy2((char
         *)&dst[idst+7*seldiff], (char *)&src7[selstart], nbytes); idst +=
         8*seldiff;
                                      }
                              }*/
      else {
        for (int ix = xstart; ix < xend; ix += 4) {
          const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));
          const Real* src1 = src0 + gptfloats;
          const Real* src2 = src1 + gptfloats;
          const Real* src3 = src2 + gptfloats;

          memcpy2(
              (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 1 * seldiff], (char*)&src1[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 2 * seldiff], (char*)&src2[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 3 * seldiff], (char*)&src3[selstart], nbytes);
          idst += 4 * seldiff;
        }
      }
    }
  }
}

inline void pack_stripes(
    const Real* const srcbase, Real* const dst, const size_t gptfloats,
    const int selstart, const int selend, const int xstart, const int ystart,
    const int zstart, const int xend, const int yend, const int zend,
    const int bx, const int by) {
  const int seldiff = selend - selstart;
  const int nbytes = seldiff * sizeof(Real);

  if ((xend - xstart) == bx) {
    for (int idst = 0, iz = zstart; iz < zend; ++iz) {
      for (int iy = ystart; iy < yend; ++iy) {
        for (int ix = xstart; ix < xend; ix += 4) {
          const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));
          const Real* src1 = src0 + gptfloats;
          const Real* src2 = src1 + gptfloats;
          const Real* src3 = src2 + gptfloats;

          memcpy2(
              (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 1 * seldiff], (char*)&src1[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 2 * seldiff], (char*)&src2[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 3 * seldiff], (char*)&src3[selstart], nbytes);
          idst += 4 * seldiff;
        }
      }
    }
  } else // bgq: 3
  {
    for (int idst = 0, iz = zstart; iz < zend; ++iz) {
      for (int iy = ystart; iy < yend; ++iy) {
        int ix = xstart;
        // for(int ix=xstart; ix<xend; ++ix)
        {
          const Real* src0 = srcbase + gptfloats * (ix + bx * (iy + by * iz));
          const Real* src1 = src0 + gptfloats;
          const Real* src2 = src1 + gptfloats;

          memcpy2(
              (char*)&dst[idst + 0 * seldiff], (char*)&src0[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 1 * seldiff], (char*)&src1[selstart], nbytes);
          memcpy2(
              (char*)&dst[idst + 2 * seldiff], (char*)&src2[selstart], nbytes);
          idst += 3 * seldiff;
        }
      }
    }
  }
}

inline void unpack(
    const Real* const pack, Real* const dstbase, const size_t gptfloats,
    const int* const selected_components, const int ncomponents,
    const int /*nsrc*/, const int dstxstart, const int dstystart,
    const int dstzstart, const int dstxend, const int dstyend,
    const int dstzend, const int xsize, const int ysize, const int /*zsize*/) {
  for (int s = 0, zd = dstzstart; zd < dstzend; ++zd)
    for (int yd = dstystart; yd < dstyend; ++yd)
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        Real* const dst =
            dstbase + gptfloats * (xd + xsize * (yd + ysize * zd));
        for (int c = 0; c < ncomponents; ++c, ++s)
          dst[selected_components[c]] = pack[s];
      }
}

inline void unpack_soa(
    const Real* const __restrict pack, Real* const __restrict dstbase,
    const int dstxstart, const int dstystart, const int dstzstart,
    const int dstxend, const int dstyend, const int dstzend, const int xsize,
    const int ysize) {
  const Real* src = pack;
  for (int zd = dstzstart; zd < dstzend; ++zd) {
    for (int yd = dstystart; yd < dstyend; ++yd) {
      Real* dst = dstbase + dstxstart + xsize * (yd + ysize * zd);
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        *dst++ = *src++;
      }
    }
  }
}

inline void unpack_soa_ncomp(
    const Real* const __restrict pack, Real* const __restrict dstbase,
    const int ncomp, const int dstxstart, const int dstystart,
    const int dstzstart, const int dstxend, const int dstyend,
    const int dstzend, const int xsize, const int ysize) {
  const Real* src = pack;
  for (int zd = dstzstart; zd < dstzend; ++zd) {
    for (int yd = dstystart; yd < dstyend; ++yd) {
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        Real* dst = dstbase + ncomp * (xd + xsize * (yd + ysize * zd));
        for (int c = 0; c < ncomp; ++c) {
          *dst++ = *src++;
        }
      }
    }
  }
}

inline void unpack1(
    const Real* const pack, Real* const dstbase, const size_t gptfloats,
    const int* const selected_components, const int ncomponents,
    const int /*nsrc*/, const int dstxstart, const int dstystart,
    const int dstzstart, const int dstxend, const int dstyend,
    const int dstzend, const int xsize, const int ysize, const int /*zsize*/) {
  int first_component = selected_components[0];
  for (int s = 0, zd = dstzstart; zd < dstzend; ++zd)
    for (int yd = dstystart; yd < dstyend; ++yd)
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        Real* const dst =
            dstbase + gptfloats * (xd + xsize * (yd + ysize * zd));
        /// for(int c=0; c<ncomponents; ++c, ++s)
        //	dst[selected_components[c]] = pack[s];
        memcpy2(
            (char*)(dst + first_component), (char*)&pack[s],
            ncomponents * sizeof(Real));
        s += ncomponents;
      }
}

inline void unpack2(
    const Real* const pack, Real* const dstbase, const size_t gptfloats,
    const int* const selected_components, const int ncomponents,
    const int /*nsrc*/, const int dstxstart, const int dstystart,
    const int dstzstart, const int dstxend, const int dstyend,
    const int dstzend, const int xsize, const int ysize, const int /*zsize*/) {
  int first_component = selected_components[0];
  const int nbytes = ncomponents * sizeof(Real);
  for (int s = 0, zd = dstzstart; zd < dstzend; ++zd)
    for (int yd = dstystart; yd < dstyend; ++yd) {
      if ((dstxend - dstxstart) % 4 != 0) {
        for (int xd = dstxstart; xd < dstxend; ++xd) {
          Real* const dst =
              dstbase + gptfloats * (xd + xsize * (yd + ysize * zd));
          // for(int c=0; c<ncomponents; ++c, ++s)
          //	dst[selected_components[c]] = pack[s];
          memcpy2((char*)(dst + first_component), (char*)&pack[s], nbytes);
          s += ncomponents;
        }
      } else {
        for (int xd = dstxstart; xd < dstxend; xd += 4) {
          Real* const dst0 = dstbase +
                             gptfloats * (xd + xsize * (yd + ysize * zd)) +
                             first_component;
          Real* const dst1 = dst0 + gptfloats;
          Real* const dst2 = dst1 + gptfloats;
          Real* const dst3 = dst2 + gptfloats;

          memcpy2((char*)dst0, (char*)&pack[s + 0 * ncomponents], nbytes);
          memcpy2((char*)dst1, (char*)&pack[s + 1 * ncomponents], nbytes);
          memcpy2((char*)dst2, (char*)&pack[s + 2 * ncomponents], nbytes);
          memcpy2((char*)dst3, (char*)&pack[s + 3 * ncomponents], nbytes);
          s += 4 * ncomponents;
        }
      }
    }
}

inline void unpack_subregion(
    const Real* const pack, Real* const dstbase, const size_t gptfloats,
    const int* const selected_components, const int ncomponents,
    const int srcxstart, const int srcystart, const int srczstart, const int LX,
    const int LY, const int dstxstart, const int dstystart, const int dstzstart,
    const int dstxend, const int dstyend, const int dstzend, const int xsize,
    const int ysize, const int /*zsize*/) {
  for (int zd = dstzstart; zd < dstzend; ++zd)
    for (int yd = dstystart; yd < dstyend; ++yd)
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        Real* const dst =
            dstbase + gptfloats * (xd + xsize * (yd + ysize * zd));
        const Real* src =
            pack + ncomponents * (xd - dstxstart + srcxstart +
                                  LX * (yd - dstystart + srcystart +
                                        LY * (zd - dstzstart + srczstart)));

        for (int c = 0; c < ncomponents; ++c)
          dst[selected_components[c]] = src[c];
      }
}

inline void unpack_subregion_soa(
    const Real* const __restrict pack, Real* const __restrict dstbase,
    const int srcxstart, const int srcystart, const int srczstart, const int LX,
    const int LY, const int dstxstart, const int dstystart, const int dstzstart,
    const int dstxend, const int dstyend, const int dstzend, const int xsize,
    const int ysize) {
  for (int zd = dstzstart; zd < dstzend; ++zd) {
    for (int yd = dstystart; yd < dstyend; ++yd) {
      Real* dst = dstbase + dstxstart + xsize * (yd + ysize * zd);
      const Real* src =
          pack + srcxstart +
          LX * (yd - dstystart + srcystart + LY * (zd - dstzstart + srczstart));
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        *dst++ = *src++;
      }
    }
  }
}

inline void unpack_subregion_soa_ncomp(
    const Real* const __restrict pack, Real* const __restrict dstbase,
    const int ncomp, const int srcxstart, const int srcystart,
    const int srczstart, const int LX, const int LY, const int dstxstart,
    const int dstystart, const int dstzstart, const int dstxend,
    const int dstyend, const int dstzend, const int xsize, const int ysize) {
  for (int zd = dstzstart; zd < dstzend; ++zd) {
    for (int yd = dstystart; yd < dstyend; ++yd) {
      Real* dst = dstbase + ncomp * (dstxstart + xsize * (yd + ysize * zd));
      for (int xd = dstxstart; xd < dstxend; ++xd) {
        const Real* src =
            pack + ncomp * (xd - dstxstart + srcxstart +
                            LX * (yd - dstystart + srcystart +
                                  LY * (zd - dstzstart + srczstart)));
        for (int c = 0; c < ncomp; ++c) {
          *dst++ = *src++;
        }
      }
    }
  }
}

#undef memcpy2
} // namespace PUPkernelsMPI

#endif /* PUPKERNELSMPI_H_XBFRP5V0 */
