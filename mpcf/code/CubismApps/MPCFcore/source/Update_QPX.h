/*
 *  Update_QPX.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/7/13.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "Update.h"

struct Update_QPX : public Update_CPP
{
	enum {
		NPOINTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_,
		JUMP = sizeof(Real) * 4
	};

public:

	Update_QPX(const Real b = 1): Update_CPP(b) {}

	template<int NFLOATS, int STEP>
	void _compute(const vector4double myb1, const vector4double myb2, Real * const src, Real * const dst) const
	{
		for(int i = 0; i < NFLOATS; i += STEP)
		{
			vec_sta(vec_madd(myb1, vec_lda(0L, src + i), vec_lda(0L, dst + i)), 0L, dst + i);
			vec_sta(vec_madd(myb2, vec_lda(JUMP, src + i), vec_lda(JUMP, dst + i)), JUMP, dst + i);
		}
	}

	void compute(const Real * const _src, Real * const dst, const int gptfloats) const
	{
		assert(gptfloats == 8);

		Real * src = const_cast<Real *>(_src);

		vector4double myb1 = vec_splats(m_b);
		vector4double myb2 = vec_perm(vec_splats(m_b), vec_splats(0), vec_gpci(01114));

		if (gptfloats == 8)
			_compute<NPOINTS * 8, 8>(myb1, myb2, src, dst);
		else
		{
			printf("ooops Update_QPX::compute: gptfloats is not quite right. aborting.\n");
			abort();
		}
	}
};
