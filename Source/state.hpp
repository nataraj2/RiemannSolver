#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "llf_flux_HLLC.hpp"

namespace Euler{

struct dynamic_state {
    std::vector<double> rho,  rhou, E;
};


class state{
    public:
        dynamic_state dynamic;
        void write_solution(const int iter);
        std::vector<double> rho_slope, rhou_slope, E_slope;
        std::vector<double> flux_rho, flux_rhou, flux_E;
        void build_u_next(dynamic_state& u_next);
        void compute_flux();
        void compute_limited_slopes_consvar();
};

void compute_slopes(const std::vector<double> &var, std::vector<double> &slope);
void boundary_conditions(dynamic_state& a_dyn_state);
void apply_bc(std::vector<double>& a_var);

void state::build_u_next(dynamic_state& u_next)
{

	compute_limited_slopes_consvar();
    compute_flux();

	for(int i=2; i<=nx-1; i++){
        u_next.rho[i]  = u_next.rho[i]  - dt*(flux_rho[i+1]-flux_rho[i])/dx;
        u_next.rhou[i] = u_next.rhou[i] - dt*(flux_rhou[i+1]-flux_rhou[i])/dx;
        u_next.E[i]    = u_next.E[i]    - dt*(flux_E[i+1]-flux_E[i])/dx;
     }

	boundary_conditions(u_next);
}


void compute_slopes(const std::vector<double> &var, std::vector<double> &slope)
{
	double deltaU_left, deltaU_right, slope1;

	for(int i=1; i<=nx; i++){
		deltaU_left =  var[i] - var[i-1];
		deltaU_right = var[i+1] - var[i];
		
		// Minmod limiter	
		slope1 = minmod(deltaU_left    , deltaU_right);
		slope[i] = slope1;	
	}	
}

void state::compute_limited_slopes_consvar()
{
	const auto &rho = dynamic.rho;
	const auto &rhou = dynamic.rhou;
	const auto &E = dynamic.E;

	rho_slope.resize(nx+2);
	rhou_slope.resize(nx+2);
	E_slope.resize(nx+2);
	
	// Compute slopes and store them 
	compute_slopes(rho , rho_slope);
	compute_slopes(rhou, rhou_slope);
	compute_slopes(E   , E_slope);	
}


void ComputeLeftAndRightStates(const std::vector<double>& a_var, const int a_n, const int a_i, const std::vector<double>& a_slope,
                                      std::vector<double>& a_state_left, std::vector<double>& a_state_right)
{
    a_state_left[a_n]  = a_var[a_i-1] + 0.5*a_slope[a_i-1];
    a_state_right[a_n] = a_var[a_i] - 0.5*a_slope[a_i];
}

void state::compute_flux()
{

	std::vector<double> QL(ND+2,0.0), QR(ND+2,0.0);

	const auto &rho = dynamic.rho;
	const auto &rhou = dynamic.rhou;
	const auto &E = dynamic.E;

	flux_rho.resize(nx+1,0.0);
	flux_rhou.resize(nx+1,0.0);
	flux_E.resize(nx+1,0.0);

	// Face loop over all faces
	for(int i=2;i<=nx; i++){

		// Compute left and right states. ie. for all the conservative variables
		ComputeLeftAndRightStates(rho,  0, i, rho_slope,  QL, QR);	
		ComputeLeftAndRightStates(rhou, 1, i, rhou_slope, QL, QR);	
		ComputeLeftAndRightStates(E   , 2, i, E_slope,    QL, QR);	

		std::array<double,3> flux = llf_flux(QL[0], QL[1], QL[2], QR[0], QR[1], QR[2]);

		flux_rho[i] = flux[0]; 
		flux_rhou[i] = flux[1]; 
		flux_E[i] = flux[2]; 
	}
}

void state::write_solution(const int iter)
{
	std::string filename;
	FILE* output_file;

	if(iter<=9){
        filename = "Solution_000" + std::to_string(iter) + ".txt";
    }
	else if(iter<=99){	
        filename = "Solution_00" + std::to_string(iter) + ".txt";
	}
	else if(iter<=999){	
        filename = "Solution_0" + std::to_string(iter) + ".txt";
	}
	else if(iter<=9999){	
        filename = "Solution_" + std::to_string(iter) + ".txt";
	}

    output_file = fopen(filename.c_str(), "w");

    for(int i=0;i<=nx+1; i++){
        fprintf(output_file, "%g %g %g %g\n",x[i], dynamic.rho[i], dynamic.rhou[i], dynamic.E[i]);
    }

    fclose(output_file);
}

void apply_bc(std::vector<double>& a_var)
{
	a_var[0]    = a_var[2];
	a_var[1]    = a_var[2];
	a_var[nx]   = a_var[nx-1];
	a_var[nx+1] = a_var[nx-1];
}

void boundary_conditions(dynamic_state& a_dyn_state)
{
	apply_bc(a_dyn_state.rho);
	apply_bc(a_dyn_state.rhou);
	apply_bc(a_dyn_state.E);
}
}
