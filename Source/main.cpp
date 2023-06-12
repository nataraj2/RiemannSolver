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

	double time = 0.0;

	int niter = t_end/dt;
	
	for(int iter=1; iter<niter; iter++){
		//dt = ComputeTimeStep(cell, CFL);

		//ComputeLimitedSlopes_ConsVar(cell);
		//ComputeLeftAndRightConsVarStatesAndComputeFlux(cell, face);
	
		ComputeLimitedSlopes_PrimVar(cell);
		ComputeLeftAndRightPrimVarStatesAndComputeFlux(cell, face);

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

