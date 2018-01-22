/*
 *  check_error.h
 *  ComputeExpansions
 *
 *  Created by Diego Rossinelli on 1/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cmath>
#include <limits>
#include <cstdio>

using namespace std;

template<typename T>
void check_error(const double tol, T ref[], T val[], const int N)
{
	static const bool verbose = false;
	int nan_ctr = 0;
	int dif_ctr = 0;

//	printf("%d\n",  N);
	for(int i=0; i<N; ++i)
	{
		assert(!std::isnan(ref[i]));
		assert(!std::isnan(val[i]));

		if (isnan(ref[i]) || isnan(val[i])) {
			nan_ctr++;
			continue;
		}
		const double err = ref[i] - val[i];
//		const double relerr = err/std::max((double)std::numeric_limits<T>::epsilon(), std::max(fabs(val[i]), fabs(ref[i])));
		const double relerr = err/std::max((double)std::numeric_limits<T>::epsilon(), (double)std::max(fabs(val[i]), fabs(ref[i])));

		if (verbose) printf("+%1.1e,", relerr);

		if (fabs(relerr) >= tol && fabs(err) >= tol) {
			dif_ctr++;
			printf("\n%d: %e %e -> %e %e\n", i, ref[i], val[i], err, relerr);
		}

		assert(fabs(relerr) < tol || fabs(err) < tol);
	}

	if (verbose) printf("\t");

	if (true) { //(verbose)
		if (nan_ctr) printf("\nnan_ctr=%d", nan_ctr);
		if (dif_ctr) printf("\ndif_ctr=%d", dif_ctr);
		if (nan_ctr || dif_ctr) printf("\n");
	}
}


template<typename T>
void check_error(const double tol, T ref, T val)
{
	check_error(tol, &ref, &val, 1);
}


// rasthofer March 2016
template<typename T>
void check_err(const double tolabs[], const double tolrel, const T ref[], const T val[], const int N, const bool verbose=true)
{
    bool bFailed = false;
  // loop all values
  for(int i=0; i<N; ++i)
  {
    // some simple tests first
    if (std::isnan(ref[i]))
    {
        bFailed = true;
       // cout << "TEST FAILED WITH NAN! REFERENCE FIELD ID: " << i << endl;
       // abort();
    }
    if (std::isinf(ref[i]))
    {
        bFailed = true;
       // cout << "TEST FAILED WITH INF! REFERENCE FIELD ID: " << i << endl;
       // abort();
    }
    if (std::isnan(val[i]))
    {
        bFailed = true;
       // cout << "TEST FAILED WITH NAN! TESTED FIELD ID: " << i << endl;
       // abort();
    }
    if (std::isinf(val[i]))
    {
        bFailed = true;
       // cout << "TEST FAILED WITH INF! TESTED FIELD ID: " << i << endl;
       // abort();
    }

    // if (verbose)
    //  printAccuracyTitle();
    // compute error
    const double error = abs((double)ref[i]-(double)val[i]);
    const double relerr = abs(error)/std::max((double)std::numeric_limits<T>::epsilon(), abs((double)ref[i]));
    // check error
    bool failed_abs = (error > tolabs[i]);
    bool failed_rel = (relerr > tolrel);
    bFailed = failed_abs || failed_rel;

    if (verbose && (failed_abs || failed_rel))
    {
      printf("\tVALUES in FIELD %d: %.12e (reference) %.12e (tested kernel)\n", i, ref[i], val[i]);
      printf("\tERROR: %e (relative: %e) (accuracy: %e abs, %e rel)\n", error, relerr, tolabs[i], tolrel);

      if (failed_abs or failed_rel)
      {
        bFailed = true;
        // cout << "\tTEST FAILED" << endl;
        // abort();
      }
    }
  }

  if (verbose && bFailed)
  {
      cout << "---> TEST FAILED!" << endl;
  }

  return;
}

template<typename T>
void check_err(const double tolabs, const double tolrel,  const T ref, const T val, const bool verbose=true)
{
  check_err(&tolabs, tolrel, &ref, &val, 1, verbose);
  return;
}
