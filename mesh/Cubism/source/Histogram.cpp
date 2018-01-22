/*
 *  Histogram.cpp
 *  
 *
 *  Created by Babak Hejazialhosseini on 3/15/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#include "Histogram.h"
#include <fstream>
#include <assert.h>
#include <math.h>
#include <iostream>

void Histogram::_setup()
{
    isroot = MPI::COMM_WORLD.Get_rank()==0;
    bInitialized = true;
}

void Histogram::notify(string sKernel, float dt)
{
    if (!bInitialized) 
        _setup();
    
    mk2t[sKernel].push_back(dt);
}

void Histogram::consolidate()
{
    assert(bInitialized);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    for(map<string,vector<float> >::iterator it=mk2t.begin(); it!=mk2t.end(); ++it)
    {
        vector<float> vSpentTime(mk2t[it->first].begin(),mk2t[it->first].end());
        
        vector<float> vSpentTime_all;
        if (isroot) vSpentTime_all.resize(vSpentTime.size()*MPI::COMM_WORLD.Get_size());
        
        MPI::COMM_WORLD.Gather(&vSpentTime.front(), vSpentTime.size(), MPI_FLOAT, &vSpentTime_all.front(), vSpentTime.size(), MPI_FLOAT, 0);
        
        if (isroot)
        {
            _print_statistcis(it->first, vSpentTime_all);
            _print2file(it->first, vSpentTime_all);
        }
    }
    
    for( map<string,vector<float> >::iterator it=mk2t.begin(); it!=mk2t.end(); ++it)
        it->second.clear();    
}

void Histogram::_print2file(string sKernel, vector<float> & buf)
{      
    char f_name_ascii[512];
    sprintf(f_name_ascii, "hist_%s", sKernel.c_str());
    
    FILE * pFile_ascii = fopen(f_name_ascii, "a");
    for(int i=0; i<buf.size(); ++i)
        fprintf(pFile_ascii, "%e\n", buf[i]);
    fclose (pFile_ascii);
}

void Histogram::_print_statistcis(string sKernel, vector<float> & buf)
{
    float sum = 0, avg = 0, std_dev = 0;
    int size = buf.size();
    
    for(int i=0; i<size; ++i)
        sum += buf[i];
    
    avg = sum/(float)size;    
    
    for(int i=0; i<size; ++i)
        std_dev += pow(buf[i]-avg,2);
    
    std_dev = sqrt(std_dev/((float)size-1));
    
    cout << sKernel << ": (Average, STD_DEV) ("<< avg << ", " << std_dev << ")" << endl;
}
