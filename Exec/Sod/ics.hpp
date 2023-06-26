#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>

#include "Riemann.H"

void Initialize(std::vector<double> &rho, std::vector<double> &rhou, std::vector<double> &E)
{

    int ND = 1;
    nx = 401;
    xmin = -10.0, xmax = 10.0;
    gamma_air = 1.4;

    double rhoL = 1.0, pL = 100000.0, uL = 0.0;
    double rhoR = 0.125, pR = 10000.0, uR = 0.0;
	
	dt = 1.25e-5;
	t_end = 0.01;

    x.resize(nx+2);

    dx = (xmax-xmin)/(nx-1);

	rho.resize(nx+2);
	rhou.resize(nx+2);
	E.resize(nx+2);

    for(int i=0;i<=nx+1;i++){
        x[i] = xmin -2.0*dx + i*dx;
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
