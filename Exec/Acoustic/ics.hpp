#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>

#include "globals.hpp"

namespace Euler{

    int ND = 1;
	const int ng = 2;
    int nx = 800;
	int sz;
    double xmin = -3.0, xmax = 3.0;
    double gamma_air = 1.4;
	double dx, dt = 1.25e-6, t_end = 0.005;

	std::vector<double> x;

void initialize(std::vector<double> &rho, std::vector<double> &rhou, std::vector<double> &E)
{

    // Number of cells in each rank;
    dx = (xmax-xmin)/nx;
    sz = nx/rank_n;

    x.resize(sz+2*ng);

    // Allocate variables with grown ghost cells
    rho.resize(sz+2*ng);
    rhou.resize(sz+2*ng);
    E.resize(sz+2*ng);

   
	// Parameters for Gaussian acoustic pulse

	double A = 20.0, S = 20.0, rho0 = 1.226, p0 = 101325.0;

	// Number of cells in each rank;
    dx = (xmax-xmin)/nx;
	sz = nx/rank_n;
	
    x.resize(sz+2*ng);

	// Allocate variables with grown ghost cells
	rho.resize(sz+2*ng);
	rhou.resize(sz+2*ng);
	E.resize(sz+2*ng);
	
	// Fill the internal region of the proc
    for(int i=ng;i<=sz+ng-1;i++){
        x[i] = xmin - ng*dx + (i+sz*rank_me+0.5)*dx;
		double pprime = A*exp(-S*x[i]*x[i]);
		double c0 = std::sqrt(gamma_air*p0/rho0);
		double uprime = 0.0;
		double rhoprime = pprime/(c0*c0);

		rho[i]  = rho0 + rhoprime;
		rhou[i] = 0.0;
		E[i]    = (p0+pprime)/(gamma_air-1.0) + 0.5*(rhou[i]*rhou[i]/rho[i]);  
	}
}
}
