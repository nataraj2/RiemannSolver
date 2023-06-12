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

    double rhoL = 1.0, pL = 1000.0, uL = 0.0;
    double rhoC = 1.0, pC = 10.0  , uC = 0.0;
    double rhoR = 1.0, pR = 100.0 , uR = 0.0;
	
	dt = 2.5e-5;

    x.resize(nx+2);
    cell.resize(nx+2);
    face.resize(nx+3);

    dx = (xmax-xmin)/(nx-1);

    for(int i=0;i<=nx+1;i++){
        auto &icell = cell[i];
        x[i] = xmin -2.0*dx + i*dx;
        if(x[i] <= 0.1){
            icell.cons_var[0] = rhoL;
            icell.cons_var[1] = rhoL*uL;
            icell.cons_var[ND+1] = pL/(gamma_air-1.0) + 0.5*rhoL*uL*uL;

        }
        else if(x[i] >=0.9){
            icell.cons_var[0] = rhoR;
            icell.cons_var[1] = rhoR*uR;
            icell.cons_var[ND+1] = pR/(gamma_air-1.0) + 0.5*rhoR*uR*uR;
        }
		else{
			icell.cons_var[0] = rhoC;
            icell.cons_var[1] = rhoC*uC;
            icell.cons_var[ND+1] = pC/(gamma_air-1.0) + 0.5*rhoC*uC*uC;
		}	
    }

    for(int i=0;i<=nx+2;i++){
        face[i].flux.resize(ND+2);
    }
}
