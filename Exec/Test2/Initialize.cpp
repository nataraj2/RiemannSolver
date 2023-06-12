#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>

#include "Riemann.H"

void Initialize()
{

    ND = 1;
    nx = 101;
    xmin = 0.0, xmax = 1.0;
    gamma_air = 1.4;

    double rhoL = 1.0, uL = -2.0, pL = 0.4;
    double rhoR = 1.0, uR =  2.0, pR = 0.4;
	
	dt = 0.0001;
	t_end = 0.15;

    x.resize(nx+2);
    cell.resize(nx+2);
    face.resize(nx+3);

    dx = (xmax-xmin)/(nx-1);

    for(int i=0;i<=nx+1;i++){
        auto &icell = cell[i];
        x[i] = xmin -2.0*dx + i*dx;
        if(x[i] <= 0.5){
            icell.cons_var[0] = rhoL;
            icell.cons_var[1] = icell.cons_var[0]*uL;
            icell.cons_var[ND+1] = pL/(gamma_air-1.0) + 0.5*rhoL*uL*uL;

        }
        else {
            icell.cons_var[0] = rhoR;
            icell.cons_var[1] = icell.cons_var[0]*uR;
            icell.cons_var[ND+1] = pR/(gamma_air-1.0) + 0.5*rhoR*uR*uR;
        }
    }

    for(int i=0;i<=nx+2;i++){
        face[i].flux.resize(ND+2);
    }
}
