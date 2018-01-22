#include <iostream>
#include <algorithm>
#include "BSplineBase.h"
#include <cassert>

class Reconstruction
{
	public:
		double DX;
		double * A; //coefficients
		
		double xmin;
		double xmax;
		
		bool OK;
		int M ; // number of coefficients

		double mean;
		int NX;
		int BC;
		double waveLength;  // Cutoff wavelength (l sub c)
		static const double BoundaryConditions[3][4];

		/*
		 Bspline constructor for reconstruction from coefficients
		*/
		 Reconstruction();

		Reconstruction(	int nx,
				double wl,
				int bc_type,
				int num_nodes,
				double xmin, double xmax);

		Reconstruction(	const double *coeff,
				int nx,
				double wl,
				int bc_type,
				int num_nodes,
				double mean, double xmin, double xmax);
		~Reconstruction();
		void Initialize( int nx,
                                double wl,
                                int bc_type,
                                int num_nodes,
                                double xmin, double xmax);

		double evaluate(double x, double * A, double mean);
		double evaluate(double x);
		double Basis (int m, double x);
		double Basis1 (int m, double x);
		double Beta (int m);
		bool Setup(int num_nodes);
		double Ratiod (int&, double &, double &);
};
