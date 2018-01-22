/*
 *  Shapes.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/24/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#ifndef SHAPES_H_W8PDLEVG
#define SHAPES_H_W8PDLEVG

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "Types.h"

//base class is a sphere
class shape
{
protected:
    double center[3], radius;
    double bbox_s[3], bbox_e[3];

    void _copy(const double c[3], const double r)
    {
        radius = r;
        for(int i=0; i<3; ++i)
        {
            center[i] = c[i];
            bbox_s[i] = c[i] - r - 1.5*static_cast<double>(Simulation_Environment::EPSILON);
            bbox_e[i] = c[i] + r + 1.5*static_cast<double>(Simulation_Environment::EPSILON);
        }
    }

public:
    shape() {}
    shape(const double c[3], const double r) { _copy(c,r); }
    shape(const shape& rhs) { _copy(rhs.center, rhs.radius); }
    shape& operator=(const shape& rhs)
    {
        if (&rhs != this) _copy(rhs.center, rhs.radius);
        return *this;
    }

    void set(const double c[3], const double r) { _copy(c,r); }

    void get_bbox(double s[3], double e[3]) const
    {
        for(int i=0; i<3; ++i)
        {
            s[i] = bbox_s[i];
            e[i] = bbox_e[i];
        }
    }

    bool inside_my_box(const double pos[3]) const
    {
        const bool bXin = pos[0]>bbox_s[0] && pos[0]<bbox_e[0];
        const bool bYin = pos[1]>bbox_s[1] && pos[1]<bbox_e[1];
        const bool bZin = pos[2]>bbox_s[2] && pos[2]<bbox_e[2];

        return bXin && bYin && bZin;
    }

    void get_center(double c[3]) const
    {
        for(int i=0; i<3; ++i)
            c[i] = center[i];
    }

    double get_rad() const
    {
        return radius;
    }

    double eval(const double pos[3]) const
    {
        double value = 0;
        // AdHoc: Change 'd < 3' to 'd < 2' to get cylinders
        #ifdef _2D_
        for (int d = 0; d < 2; d++) {
        #else
        for (int d = 0; d < 3; d++) {
        #endif
          if (Simulation_Environment::BC_PERIODIC[d])
          {
            const double tmp = std::min(std::abs(pos[d] - (center[d]+static_cast<double>(Simulation_Environment::extents[d]))),std::abs(pos[d] - (center[d]-static_cast<double>(Simulation_Environment::extents[d]))));
            value += std::pow(std::min(abs(pos[d]-center[d]), tmp), 2.0);
          }
          else
            value += std::pow(pos[d]-center[d], 2.0);
        }

        return std::sqrt(value) - radius;
    }

    //Every other derived shape should implement this method.
    bool rejection_check(shape * this_shape, const double start[3], const double end[3]) const
    {
        double s[3], e[3];
        this->get_bbox(s,e);

        //this rule checks that the buble is inside the bounding box
        const bool bOut = s[0]<start[0] || s[1]<start[1] || s[2]<start[2] ||
        e[0]>end[0] || e[1]>end[1] || e[2]>end[2];

        if (bOut)
            return true;

        if(this!=this_shape)
        {
            double this_s[3], this_e[3];

            this_shape->get_bbox(this_s,this_e);

            const double overlap_start[3] =
            {
                std::max(this_s[0], s[0]),
                std::max(this_s[1], s[1]),
                std::max(this_s[2], s[2])
            };

            const double overlap_end[3] =
            {
                std::min(this_e[0], e[0]),
                std::min(this_e[1], e[1]),
                std::min(this_e[2], e[2])
            };

            const bool bOverlap = overlap_end[0] > overlap_start[0] && overlap_end[1] > overlap_start[1] && overlap_end[2] > overlap_start[2];

            if (bOverlap)
                return true;
        }

        return false;
    }

    static std::vector<shape> make_many(const double h, string filename, const bool verbose=true)
    {
        std::vector<shape> v_shapes;
        std::ifstream f_read_cloud(filename.c_str());

        if (!f_read_cloud.good())
        {
            std::cout << "Watchout! cant read the file " << filename << ". Aborting now..." << std::endl;
            abort();
        }

        while (true) {
            if (!f_read_cloud.good()) abort();

            int idx;
            double c[3], rad;

            // more bubbles there?
            f_read_cloud >> idx >> c[0] >> c[1] >> c[2] >> rad;
            if (f_read_cloud.tellg() == -1) break;
            else
            {
                v_shapes.push_back( shape(c,rad) );
                if (verbose) std::cout << "shape " << idx << " " <<  c[0] << " " << c[1] << " " << c[2] << " " << rad << std::endl;
            }
        }

        f_read_cloud.close();
        return v_shapes;
    }

};


template<class Tshape=shape>
class Seed
{
    std::vector<Tshape> v_shapes;

public:
    Seed() {}
    Seed(const Seed& rhs) : v_shapes(rhs.v_shapes) {}
    Seed& operator=(const Seed& rhs)
    {
        if (this != &rhs)
            v_shapes = rhs.v_shapes;
        return *this;
    }

    void make_shapes(string filename, const double h, const bool verbose=true)
    {
        v_shapes = Tshape::make_many(h, filename, verbose);
	    if (verbose) std::cout << "number of shapes are " << v_shapes.size() << std::endl;
    }

    /*------------------------------------------------------------------------*
     * identify bubbles that have to be known by each subdomain (block/rank)
     *                    revised with respect to periodic boundary conditions
     *                                                     rasthofer June 2015
     *------------------------------------------------------------------------*/
    Seed retain_shapes(const double mystart[3], const double extent[3]) const
    {
        assert(v_shapes.size() > 0);
        //v_shapes.clear();


        // get bounding box of considered subdomain
        const double myend[3] = {
            mystart[0] + extent[0],
            mystart[1] + extent[1],
            mystart[2] + extent[2]};

#ifdef _KEEPALL_
        // issue: for initial pressure field according to Tiwari, the region of influnce of an individual
        // bubble is infinity due to the application of tanh
        // there are basically 2 options to fix this
        // I> keep all bubbles on all ranks/blocks
        // 2> clip tanh and extend the bounding box of an individual bubble accordingly

        Seed<Tshape> retval;
        for(int i=0; i<v_shapes.size(); ++i)
            retval.v_shapes.push_back(v_shapes[i]);
#else

        Seed<Tshape> retval;

        // loop all bubbles
        for(int i=0; i<v_shapes.size(); ++i)
        {
            Tshape curr_shape = v_shapes[i];

            // get bounding box of current bubble
            double s[3],e[3];
            curr_shape.get_bbox(s,e);

            // check for overlap of bounding box of current bubble with bounding box
            // of evaluated subdomain
            std::vector<double> range(3,0.0);
            for (int d = 0; d < 3; d++)
                range[d] = std::min(myend[d], e[d]) - std::max(mystart[d], s[d]);

            bool bOverlap = (range[0] > 0.0) && (range[1] > 0.0) && (range[2] > 0.0);

            // if there is no overlap check also for pontential periodic boundary conditions
            // that is, shift the bounding box of current bubble with respect to the periodic
            // directions and check for overlaps once again
            if (not bOverlap)
            {
                if (Simulation_Environment::BC_PERIODIC [0] or
                        Simulation_Environment::BC_PERIODIC [1] or
                        Simulation_Environment::BC_PERIODIC [2] )
                {
                    // loop all dimensions
                    for (int d = 0; d < 3; d++)
                    {
                        if (Simulation_Environment::BC_PERIODIC [d])
                        {
                            // initialize shifted bounding box
                            double try_s = s[d];
                            double try_e = e[d];
                            // check for bounding boxes exceeding the domain
                            // and shift the box with respect to the periodic domain
                            // such that this part enters the domain from the other side
                            if (s[d] < 0.0)
                            {
                                try_s += static_cast<double>(Simulation_Environment::extents[d]);
                                try_e += static_cast<double>(Simulation_Environment::extents[d]);
                            }
                            else {
                                if (e[d] > Simulation_Environment::extents[d])
                                {
                                    try_s -= static_cast<double>(Simulation_Environment::extents[d]);
                                    try_e -= static_cast<double>(Simulation_Environment::extents[d]);
                                }
                            }

                            // check for overlap of shifted bounding box of current bubble
                            // with bounding box of evaluated subdomain
                            range[d] = std::min(myend[d], try_e) - std::max(mystart[d], try_s);

                            bOverlap = (range[0] > 0.0) && (range[1] > 0.0) && (range[2] > 0.0);

                            if (bOverlap)
                                break;

                            // revert shift if value is out of the domain
                            if (range[d]<0.0)
                                range[d] = std::min(myend[d], e[d]) - std::max(mystart[d], s[d]);
                        }
                    }
                }
            }

            // store bubble
            if (bOverlap) retval.v_shapes.push_back(curr_shape);
        }
#endif
        return retval;
    }

    std::vector<Tshape> get_shapes() const { return v_shapes; }
    std::vector<Tshape>& get_shapes() { return this->v_shapes; }
    int get_shapes_size() const { return v_shapes.size(); }
    void set_shapes(const std::vector<Tshape>& v) { v_shapes = v; }
};


template<class Tshape>
inline double distance(const std::vector<Tshape>& v_shapes, const double pos[3])
{
    double d = (double)HUGE_VAL;

    for( int i=0; i<v_shapes.size(); ++i)
    {
        const double newdistance = v_shapes[i].eval(pos);
        d = std::min(d, newdistance);
    }

    return d;
}


template<class Tshape>
inline double eval(const std::vector<Tshape>& v_shapes, const double pos[3])
{
    // determine alpha2 over shape interface

    const double d = distance(v_shapes, pos);

    const double eps = Simulation_Environment::EPSILON;
    const double alpha = M_PI*std::min(1., std::max(0., (double)(d + eps)/(2 * eps)));

    return 0.5 + 0.5 * std::cos(alpha);
}

template<class Tshape>
inline double eval_shifted(const std::vector<Tshape>& v_shapes, const double pos[3])
{
    const double d = distance(v_shapes, pos);

    const double eps = Simulation_Environment::EPSILON;
    const double alpha = M_PI*std::min(1., std::max(0., (double)(d + eps)/(2 * eps)));//this 0 shift for pressure is very important.

    return 0.5 + 0.5 * std::cos(alpha);
}

template<class Tshape>
inline double eval(const std::vector<Tshape>& v_shapes, const double pos[3], const double epsh)
{
    // determine alpha2 over shape interface

    const double d = distance(v_shapes, pos);
    const double alpha = M_PI*std::min(1., std::max(0., (double)(d + epsh)/(2 * epsh)));

    return 0.5 + 0.5 * std::cos(alpha);
}

#endif /* SHAPES_H_W8PDLEVG */
