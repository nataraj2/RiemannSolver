#include <iostream>
#include <string>
#include <vector>

#include "ics.hpp"
#include "state.hpp"

int nx, ny;
double xmin, xmax;
double gamma_air;
double dx, dt, t_end;

std::vector<double> x;

int main(){

	int freq = 10;
	double CFL = 0.25;

	state state_;

	std::vector<double> rho, rhou, E;

	Initialize(rho, rhou, E);
	state_.dynamic.rho = rho;
	state_.dynamic.rhou = rhou;
	state_.dynamic.E = E;
	state_.WriteSolution(0);
	
	double time = 0.0;
	int niter = t_end/dt;
	
	// Time advance
	for(int iter=1; iter<=niter; iter++){

		state_.build_u_next(state_.dynamic);

		if(iter%freq==0){
			state_.WriteSolution(iter/freq);
		}
		time = time + dt;
		std::cout << "Iteration is " << iter << " dt = " << dt << " " << "time = " << time << "\n";
	}
	return 0;

}

