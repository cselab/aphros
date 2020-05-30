/**
 * @file interface.c
 * @authors Simone Bnà, Sandro Manservisi, Ruben Scardovelli,
 *          Philip Yecko and Stephane Zaleski
 * @date  12 November 2015
 * @brief FORTRAN to C interface for the functions vofi_Get_fh
 *        vofi_Get_cc.
 */

#include "vofi.h"
#include "vofi_stddecl.h"

/* ------------------------------------------------------------------- *
 * DESCRIPTION:                                                        *
 * FORTRAN to C interface for the function vofi_Get_fh                      *
 * INPUT and OUTPUT: see vofi_Get_fh                                        *
 * ------------------------------------------------------------------- */

vofi_real EXPORT(vofi_get_fh)(
    integrand impl_func, vofi_creal x0[], vofi_creal* H0, vofi_cint* Ndim0,
    vofi_cint* iX0) {
  vofi_creal h0 = *H0;
  vofi_cint ndim0 = *Ndim0, ix0 = *iX0;
  vofi_real Fh;

  Fh = vofi_Get_fh(impl_func, x0, h0, ndim0, ix0);

  return Fh;
}

/* ------------------------------------------------------------------- *
 * DESCRIPTION:                                                        *
 * FORTRAN to C interface for the function vofi_Get_cc                      *
 * INPUT and OUTPUT: see vofi_Get_cc                                        *
 * ------------------------------------------------------------------- */

vofi_real EXPORT(vofi_get_cc)(
    integrand impl_func, vofi_creal x0[], vofi_creal* H0, vofi_creal* Fh,
    vofi_cint* Ndim0) {
  vofi_creal h0 = *H0, fh = *Fh;
  vofi_cint ndim0 = *Ndim0;
  vofi_real CC;

  CC = vofi_Get_cc(impl_func, x0, h0, fh, ndim0);

  return CC;
}
