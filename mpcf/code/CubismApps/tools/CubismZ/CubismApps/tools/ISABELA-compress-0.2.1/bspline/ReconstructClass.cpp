#include <iostream>
#include <algorithm>
#include "BSplineBase.h"
#include <cassert>
#include "ReconstructClass.h"
#include <cmath>

//using namespace std;
Reconstruction::Reconstruction(){}

Reconstruction::Reconstruction(	int nx,
				double wl,
				int bc_type,
				int num_nodes,
				double xmin, double xmax)
{
	A = NULL;
	M = num_nodes;
	NX = nx;
	this->xmin = xmin;
	this->xmax = xmax;
	waveLength = wl;

	BC = bc_type;

	OK = Setup(num_nodes);
}

Reconstruction::Reconstruction(const double *coeff,
                                      int nx,                                      
                                      double wl,
                                      int bc_type,
                                      int num_nodes,
									  double mean, double xmin, double xmax) 	
	{
		A =  new double[num_nodes];

		OK=true;
		M = num_nodes;
		NX = nx;
		this->xmin = xmin;
		this->xmax = xmax;
		waveLength = wl;

		//DX = (xmax - xmin) / M;

		this->mean= mean;
		BC = bc_type;

		for(int i =0 ; i < num_nodes; i++)
		{
			A[i] = coeff[i];
		}

		OK = Setup(num_nodes);
	}

void Reconstruction::Initialize( int nx,
                                double wl,
                                int bc_type,
                                int num_nodes,
                                double xmin, double xmax)
{
        A = NULL;
        M = num_nodes;
        NX = nx;
        this->xmin = xmin;
        this->xmax = xmax;
        waveLength = wl;

        BC = bc_type;

        OK = Setup(num_nodes);
}

	
Reconstruction::~Reconstruction()
{
	if(A!=NULL)
		free(A);
}

bool Reconstruction::Setup(int num_nodes)
{

    int ni = 9; // Number of node intervals (NX - 1)
    double deltax;

    if (num_nodes >= 2) {
        // We've been told explicitly the number of nodes to use.
        ni = num_nodes - 1;
        if (waveLength == 0) {
            waveLength = 1.0;
        }
    } else if (waveLength == 0) {
        // Turn off frequency constraint and just set two node intervals per
        // data point.
        ni = NX * 2;
        waveLength = 1;
    } else if (waveLength > xmax - xmin) {
        return (false);
    } else {
        // Minimum acceptable number of node intervals per cutoff wavelength.
        static const double fmin = 2.0;

        double ratiof; // Nodes per wavelength for current deltax
        double ratiod; // Points per node interval

        // Increase the number of node intervals until we reach the minimum
        // number of intervals per cutoff wavelength, but only as long as 
        // we can maintain at least one point per interval.
        do {
            if (Ratiod(++ni, deltax, ratiof) < 1.0)
                return false;
        } while (ratiof < fmin);

        // Now increase the number of intervals until we have at least 4
        // intervals per cutoff wavelength, but only as long as we can
        // maintain at least 2 points per node interval.  There's also no
        // point to increasing the number of intervals if we already have
        // 15 or more nodes per cutoff wavelength.
        // 
        do {
            if ((ratiod = Ratiod(++ni, deltax, ratiof)) < 1.0 || ratiof > 15.0) {
                --ni;
                break;
            }
        } while (ratiof< 4 || ratiod> 2.0);
    }

    // Store the calculations in our state
    M = ni;
    DX = (xmax - xmin) / ni;
	return(true);
}

//////////////////////////////////////////////////////////////////////
double Reconstruction::evaluate(double x, double *A, double mean) 
{
    double y = 0;
    if (OK) {
        int n = (int)((x - xmin)/DX);
        for (int i = std::max(0, n-1); i <= std::min(M, n+2); ++i) {
            y += A[i] * Basis(i, x);
		}
		
        y += mean;
    }

    return y;
}

//////////////////////////////////////////////////////////////////////
double Reconstruction::evaluate(double x) 
{
    double y = 0;

    if (OK) {
        int n = (int)((x - xmin)/DX);
        for (int i = std::max(0, n-1); i <= std::min(M, n+2); ++i) {
            y += A[i] * Basis(i, x);
		}
		
        y += mean;
    }

    return y;
}

double Reconstruction::Basis1(int m,double x)
{
    double y = 0;
    double xm = xmin + (m * DX);
    double z = std::abs((double)(x - xm) / (double)DX);
    if (z < 2.0) {
        z = 2 - z;
        y = 0.25 * (z*z*z);
        z -= 1.0;
        if (z > 0)
            y -= (z*z*z);
    }
   return y;
}

double Reconstruction::Basis(int m,double x)
{
    double y = 0;
    double xm = xmin + (m * DX);
    double z = std::abs((double)(x - xm) / (double)DX);
    if (z < 2.0) {
        z = 2 - z;
        y = 0.25 * (z*z*z);
        z -= 1.0;
        if (z > 0)
            y -= (z*z*z);
    }
    // Boundary conditions, if any, are an additional addend.
    if (m == 0 || m == 1)
        y += Beta(m) * Basis(-1, x);
    else if (m == M-1 || m == M)
        y += Beta(m) * Basis(M+1, x);

    return y;
}

inline double Reconstruction::Beta(int m)
{
    if (m > 1 && m < M-1)
        return 0.0;
    if (m >= M-1)
        m -= M-3;
    assert(0 <= BC && BC <= 2);
    assert(0 <= m && m <= 3);
    return BoundaryConditions[BC][m];
}

inline double Reconstruction::Ratiod(int &ni,double &deltax, double &ratiof)
{
    deltax = (xmax - xmin) / ni;
    ratiof = waveLength / deltax;
    double ratiod = (double) NX / (double) (ni + 1);
    return ratiod;
}

//////////////////////////////////////////////////////////////////////

// This array contains the beta parameter for the boundary condition
// constraints.  The boundary condition type--0, 1, or 2--is the first
// index into the array, followed by the index of the endpoints.  See the
// Beta() method.
const double Reconstruction::BoundaryConditions[3][4] =
    {
            //  0       1       M-1     M
                {
                        -4,
                        -1,
                        -1,
                        -4 },
                {
                        0,
                        1,
                        1,
                        0 },
                {
                        2,
                        -1,
                        -1,
                        2 } };

//////////////////////////////////////////////////////////////////////
