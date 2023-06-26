#include <iostream>
#include <string>
#include <vector>

#include "ics.hpp"
#include "state.hpp"

using namespace Euler;

int main(){

	int freq = 10;
	double CFL = 0.25;

	state state_;

	std::vector<double> rho, rhou, E;

	initialize(rho, rhou, E);
	state_.dynamic.rho = rho;
	state_.dynamic.rhou = rhou;
	state_.dynamic.E = E;
	state_.write_solution(0);
	
	double time = 0.0;
	int niter = t_end/dt;
	
	// Time advance
	for(int iter=1; iter<=niter; iter++){

		state_.build_u_next(state_.dynamic);

		if(iter%freq==0){
			state_.write_solution(iter/freq);
		}
		time = time + dt;
		std::cout << "Iteration is " << iter << " dt = " << dt << " " << "time = " << time << "\n";
	}
	return 0;

}

