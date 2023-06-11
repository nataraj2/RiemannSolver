#include <iostream>
#include <string>
#include <vector>

#include "Riemann.H"

int main(){

	int freq = 10;
	double CFL = 0.25;

	Initialize();
	WriteSolution(cell, 0);
	
	// Time advance

	double t_end = 0.02, time = 0.0;
	

	for(int iter=1; iter<=400; iter++){
		double dt = 2.5e-5;ComputeTimeStep(cell, CFL);
		ComputeLimitedSlopes_ConsVar(cell);
		ComputeLeftAndRightStatesAndComputeFlux(cell, face);
		std::cout << "Iteration is " << iter << " dt = " << dt << " " << "time = " << time << "\n";
		for(int i=2; i<=nx-1; i++){
			for(int n=0; n<=ND+1; n++){
				cell[i].cons_var[n] = cell[i].cons_var[n] - dt*(face[i+1].flux[n]-face[i].flux[n])/dx;
			}
		}
		BoundaryConditions(cell);
		if(iter%freq==0){
			WriteSolution(cell, iter/freq);
		}
		time = time + dt;
	}
	return 0;

}

