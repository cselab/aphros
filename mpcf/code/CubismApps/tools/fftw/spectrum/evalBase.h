/* File:   evalBase.h */
/* Date:   July 2015 */
/* Author: Ursula Rasthofer */
/* Tag:    base class for evaluations in spectral space */
/* Copyright 2015 ETH Zurich. All Rights Reserved. */

#ifndef EVALBASE_H
#define EVALBASE_H

#include "types.h"

class evalBase
{

public:
    virtual ~evalBase() {}
    virtual void compute() = 0;
};

#endif
