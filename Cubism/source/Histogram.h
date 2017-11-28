/*
 *  Histogram.h
 *  
 *
 *  Created by Babak Hejazialhosseini on 3/15/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <mpi.h>
#include <vector>
#include <map>
#include <string>

using namespace std;

class Histogram
{
    map<string,vector<float> > mk2t;
    bool isroot;
    bool bInitialized;
    int reportID;
    
    void _print2file(string sKernel, vector<float> & buf);
    void _print_statistcis(string sKernel, vector<float> & buf);
    void _setup();
    
public:
    Histogram(int a=0): mk2t() {}
    
    void notify(string sKernel, float dt);
    void consolidate();
};
