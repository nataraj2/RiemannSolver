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
    double xmin = -10.0, xmax = 10.0;
    double gamma_air = 1.4;
	double dx, dt = 1.25e-5/2.0, t_end = 0.01+1e-10;

	std::vector<double> x;

void initialize(std::vector<double> &rho, std::vector<double> &rhou, std::vector<double> &E)
{
    double rhoL = 1.0, pL = 100000.0, uL = 0.0;
    double rhoR = 0.125, pR = 10000.0, uR = 0.0;

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
        if(x[i] <= 0.0){
            rho[i] = rhoL;
            rhou[i] = rhoL*uL;
            E[i] = pL/(gamma_air-1.0) + 0.5*rhoL*uL*uL;

        }
        else if (x[i]>=0.0){
            rho[i]  = rhoR;
            rhou[i] = rhoR*uR;
            E[i] = pR/(gamma_air-1.0) + 0.5*rhoR*uR*uR;
        }
    }
}
}
