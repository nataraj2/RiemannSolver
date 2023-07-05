#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "llf_flux_HLLC.hpp"

#include "mpi.h"

namespace Euler{

struct msg_t {
    msg_t()=default;
    msg_t(const std::array<std::array<double,3>,ng>& arr)
	{
		rho_msg.resize(ng,0.0);
		rhou_msg.resize(ng,0.0);
		E_msg.resize(ng,0.0);

		for(int i=0; i<ng; i++){
			rho_msg[i]  = arr[i][0];
			rhou_msg[i] = arr[i][1];
			E_msg[i]    = arr[i][2];
		}	   
	}

    void write_to_array(std::array<std::array<double,3>,ng>& arr) {
		for (int i = 0; i < ng; ++i){
        	arr[i][0] = rho_msg[i];
        	arr[i][1] = rhou_msg[i];
        	arr[i][2] = E_msg[i];
    	}
	}

	std::vector<double> rho_msg, rhou_msg, E_msg;
	std::array<std::array<double, 3>,ng> send_buffer, receive_buffer;
};


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
        void lift_bdry(dynamic_state& u_next);
        void compute_flux();
        void compute_limited_slopes_consvar();
};

void compute_slopes(const std::vector<double> &var, std::vector<double> &slope);
void boundary_conditions(dynamic_state& a_dyn_state);
void apply_bc(std::vector<double>& a_var);

/*inline state::msg_t state::make_message(uint id_me, uint id_to) {

    return id_me > id_to ?
        msg_t(dynamic.rho[0], dynamic.rhou[0], dynamic.E[0]) :
        msg_t(dynamic.rho[sz-1], dynamic.rhou[sz-1], dynamic.E[sz-1]);
}*/

void state::build_u_next(dynamic_state& u_next)
{

	compute_limited_slopes_consvar();
    compute_flux();

	for(int i=2*ng; i<=sz-1; i++){
        u_next.rho[i]  = u_next.rho[i]  - dt*(flux_rho[i+1]-flux_rho[i])/dx;
        u_next.rhou[i] = u_next.rhou[i] - dt*(flux_rhou[i+1]-flux_rhou[i])/dx;
        u_next.E[i]    = u_next.E[i]    - dt*(flux_E[i+1]-flux_E[i])/dx;
     }

	boundary_conditions(u_next);
}

void state::lift_bdry(dynamic_state& u_next)
{

	for(int i=ng; i<=2*ng-1; i++){
        u_next.rho[i]  = u_next.rho[i]  - dt*(flux_rho[i+1]-flux_rho[i])/dx;
        u_next.rhou[i] = u_next.rhou[i] - dt*(flux_rhou[i+1]-flux_rhou[i])/dx;
        u_next.E[i]    = u_next.E[i]    - dt*(flux_E[i+1]-flux_E[i])/dx;
     }

	for(int i=sz; i<=sz+ng-1; i++){
        u_next.rho[i]  = u_next.rho[i]  - dt*(flux_rho[i+1]-flux_rho[i])/dx;
        u_next.rhou[i] = u_next.rhou[i] - dt*(flux_rhou[i+1]-flux_rhou[i])/dx;
        u_next.E[i]    = u_next.E[i]    - dt*(flux_E[i+1]-flux_E[i])/dx;
     }
}

void compute_slopes(const std::vector<double> &var, std::vector<double> &slope, const int i)
{
	double deltaU_left, deltaU_right, slope1;

	deltaU_left =  var[i] - var[i-1];
	deltaU_right = var[i+1] - var[i];
		
	// Minmod limiter	
	slope1 = minmod(deltaU_left    , deltaU_right);
	slope[i] = slope1;	
}

void state::compute_limited_slopes_consvar()
{
	const auto &rho = dynamic.rho;
	const auto &rhou = dynamic.rhou;
	const auto &E = dynamic.E;

	rho_slope.resize(sz+2*ng);
	rhou_slope.resize(sz+2*ng);
	E_slope.resize(sz+2*ng);
	
	// Compute slopes and store them 
	for(int i=ng-1; i<=sz+ng; i++){
		compute_slopes(rho , rho_slope, i);
		compute_slopes(rhou, rhou_slope, i);
		compute_slopes(E   , E_slope, i);	
	}	
}


void compute_left_and_right_states(const std::vector<double>& a_var, const int a_n, const int a_i, const std::vector<double>& a_slope,
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

	flux_rho.resize(sz+2*ng+1,0.0);
	flux_rhou.resize(sz+2*ng+1,0.0);
	flux_E.resize(sz+2*ng+1,0.0);

	// Face loop over all faces
	for(int i=ng;i<=sz+ng; i++){

		// Compute left and right states. ie. for all the conservative variables
		compute_left_and_right_states(rho,  0, i, rho_slope,  QL, QR);	
		compute_left_and_right_states(rhou, 1, i, rhou_slope, QL, QR);	
		compute_left_and_right_states(E   , 2, i, E_slope,    QL, QR);	

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

    output_file = fopen(filename.c_str(), "a");

	for(int i = 0; i<rank_n; i++){
		if(rank_me == i){	
    		for(int i=ng;i<=sz+ng-1; i++){
        		fprintf(output_file, "%g %g %g %g\n",x[i], dynamic.rho[i], dynamic.rhou[i], dynamic.E[i]);
    		}
			fflush(output_file);
		}	
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

    fclose(output_file);
}

void apply_bc(std::vector<double>& a_var)
{
	for(int i=0;i<=ng-1;i++){
		a_var[i]    = a_var[ng];
	}	
	for(int i=sz+ng;i<=sz+2*ng-1;i++){
		a_var[i]   = a_var[sz+ng-1];
	}
}

void boundary_conditions(dynamic_state& a_dyn_state)
{
	apply_bc(a_dyn_state.rho);
	apply_bc(a_dyn_state.rhou);
	apply_bc(a_dyn_state.E);
}
}
