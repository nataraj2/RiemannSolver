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
		rho.resize(ng,0.0);
		rhou.resize(ng,0.0);
		E.resize(ng,0.0);

		for(int i=0; i<=ng-1; i++){
			rho[i]  = arr[i][0];
			rhou[i] = arr[i][1];
			E[i]    = arr[i][2];
		}	   
	}

    inline void write_to_array(std::array<std::array<double,3>,ng>& arr) {
		for (int i = 0; i<=ng-1; ++i){
        	arr[i][0] = rho[i];
        	arr[i][1] = rhou[i];
        	arr[i][2] = E[i];
    	}
	}

	std::vector<double> rho, rhou, E;
	std::array<std::array<double, 3>,ng> send_buffer, receive_buffer;
};

struct bdry_t {
	bdry_t(std::array<std::array<double,3>,ng>& arr){	
		flux_diff_rho.resize(ng,0.0);	
		flux_diff_rhou.resize(ng,0.0);	
		flux_diff_E.resize(ng,0.0);

		for(int i=0; i<=ng-1; i++){
			flux_diff_rho[i]  = arr[i][0];
			flux_diff_rhou[i] = arr[i][1];
			flux_diff_E[i]    = arr[i][2];
		}
	}		
	std::vector<double> flux_diff_rho, flux_diff_rhou, flux_diff_E;		
};

inline bdry_t operator * (double a, const bdry_t& bdry) {

	std::array<std::array<double,3>,ng> arr;
	for(int i=0; i<=ng-1; i++){
		arr[i][0] = a*bdry.flux_diff_rho[i];	
		arr[i][1] = a*bdry.flux_diff_rhou[i];	
		arr[i][2] = a*bdry.flux_diff_E[i];	
	}		
    return bdry_t(arr);
}

struct dynamic_state {
    std::vector<double> rho,  rhou, E;
};

enum Side : unsigned char {left = 0, right = 1};

class state{
	using msg_t = Euler::msg_t;
	using bdry_t = Euler::bdry_t;
    public:
        dynamic_state dynamic;
        std::vector<double> rho_slope, rhou_slope, E_slope;
        std::vector<double> flux_rho, flux_rhou, flux_E;
		std::vector<Side> physical_boundaries;
		void allocate_arrays();
		msg_t make_message(uint id_from , uint id_to);
        void build_u_next(dynamic_state& u_next);
		state::bdry_t process_message(uint id_me, uint id_from, const msg_t& msg );
        void lift_bdry(uint id_me, uint id_from, dynamic_state& u_next, const bdry_t& bdry);
        void compute_flux(int i);
        void compute_limited_slopes_consvar(int i);
		void boundary_conditions(dynamic_state& a_dyn_state);
		void apply_bc(std::vector<double>& a_var);
        void write_solution(const int iter);
};

void compute_slopes(const std::vector<double> &var, std::vector<double> &slope);
void boundary_conditions(dynamic_state& a_dyn_state);
void apply_bc(std::vector<double>& a_var);

inline void state::allocate_arrays()
{
	rho_slope.resize(sz+2*ng);
	rhou_slope.resize(sz+2*ng);
	E_slope.resize(sz+2*ng);

	flux_rho.resize(sz+2*ng+1,0.0);
	flux_rhou.resize(sz+2*ng+1,0.0);
	flux_E.resize(sz+2*ng+1,0.0);
}	

inline void state::build_u_next(dynamic_state& u_next)
{
	// Do all internal cells that do not require communicated data

	// Do boundary conditions
	boundary_conditions(u_next);
	
	// Rank 0 and n-1 can do extra cells during this stage
	int istart, iend;
	if(physical_boundaries.size()>0){
		for(auto bdry: physical_boundaries){
			if(bdry==left){
				istart = 2;
				iend = sz-1;
			}
			else if(bdry==right){
				istart = 2*ng;
				iend = sz+ng-1;
			}
			else{
				std::cout << "Should not reach here in build_u_next. Exiting...." << "\n";
				exit(1);
			}
		}
	}
	else{
		istart = 2*ng;
		iend = sz-1;
	}

	// Compute all internal slopes
	for(int i=istart-1; i<=iend+1; i++){
		compute_limited_slopes_consvar(i);
	}

	// Compute all internal fluxes
	for(int i=istart;i<=iend+1; i++){
    	compute_flux(i);
	}

	// Compute all internal updates
	for(int i=istart; i<=iend; i++){
        u_next.rho[i]  = u_next.rho[i]  - dt*(flux_rho[i+1]-flux_rho[i])/dx;
        u_next.rhou[i] = u_next.rhou[i] - dt*(flux_rhou[i+1]-flux_rhou[i])/dx;
        u_next.E[i]    = u_next.E[i]    - dt*(flux_E[i+1]-flux_E[i])/dx;
     }
}

inline state::msg_t state::make_message(uint id_me, uint id_to) {

	std::array<std::array<double,3>,ng> arr;

    if(id_me > id_to){
		for(int i=0; i<=ng-1; i++){
			arr[i][0] = dynamic.rho[ng+i];				
			arr[i][1] = dynamic.rhou[ng+i];				
			arr[i][2] = dynamic.E[ng+i];
		}
	}
	else if (id_me < id_to) {
		for(int i=0; i<=ng-1; i++){
			arr[i][0] = dynamic.rho[sz+i];				
			arr[i][1] = dynamic.rhou[sz+i];				
			arr[i][2] = dynamic.E[sz+i];
		}
	}
	else{
		std::cout << "Should not reach here inside make message. Exiting ... " << "\n";
		exit(1);
	}
	return msg_t(arr);
}

inline state::bdry_t state::process_message(uint id_me, uint id_from, const msg_t& msg )
{	
	std::array<std::array<double,3>,ng> arr;

	auto &rho = dynamic.rho;
	auto &rhou = dynamic.rhou;
	auto &E = dynamic.E;
	
	if(id_me > id_from){		
		for(int i=0; i<=ng-1; i++){
			rho[i] = msg.rho[i];
			rhou[i] = msg.rhou[i];
			E[i] = msg.E[i];
		}
		for(int i=ng-1; i<=2*ng-2; i++){
	    	compute_limited_slopes_consvar(i);
		}

		for(int i=ng; i<=2*ng-1; i++){
    		compute_flux(i);
		}
		for(int i=0;i<=ng-1;i++){
			arr[i][0] = (flux_rho[ng+1+i]-flux_rho[ng+i]);
			arr[i][1] = (flux_rhou[ng+1+i]-flux_rhou[ng+i]);
			arr[i][2] = (flux_E[ng+1+i]-flux_E[ng+i]);
		}
	}
	else if(id_me < id_from){
		for(int i=0; i<=ng-1; i++){
            rho[sz+ng+i] = msg.rho[i];
            rhou[sz+ng+i] = msg.rhou[i];
            E[sz+ng+i] = msg.E[i];
        }

		for(int i=sz+1; i<=sz+ng; i++){
    		compute_limited_slopes_consvar(i);
		}

		for(int i=sz+1; i<=sz+ng; i++){
    		compute_flux(i);
		}

		for(int i=0;i<=ng-1;i++){
			arr[i][0] = (flux_rho[sz+1+i]-flux_rho[sz+i]);
			arr[i][1] = (flux_rhou[sz+1+i]-flux_rhou[sz+i]);
			arr[i][2] = (flux_E[sz+1+i]-flux_E[sz+i]);
		}
	}
	else{
		std::cout << "Should not reach here in process message. Exiting...." << "\n";
		exit(1);
	}

	return bdry_t(arr);
}


inline void state::lift_bdry(uint id_me, uint id_from, dynamic_state& u_next, const bdry_t& bdry)
{

	// Do all processor boundary cells - i.e. ng cells will be updated on each side
	if(id_me > id_from){
		for(int i=ng; i<=2*ng-1; i++){
        	u_next.rho[i]  = u_next.rho[i]  - bdry.flux_diff_rho[i-ng]/dx;
        	u_next.rhou[i] = u_next.rhou[i] - bdry.flux_diff_rhou[i-ng]/dx;
        	u_next.E[i]    = u_next.E[i]    - bdry.flux_diff_E[i-ng]/dx;
     	}
	}
	else if(id_me < id_from){
		for(int i=sz; i<=sz+ng-1; i++){
        	u_next.rho[i]  = u_next.rho[i]  - bdry.flux_diff_rho[i-sz]/dx;
        	u_next.rhou[i] = u_next.rhou[i] - bdry.flux_diff_rhou[i-sz]/dx;
        	u_next.E[i]    = u_next.E[i]    - bdry.flux_diff_E[i-sz]/dx;
     	}
	}
	else{
		std::cout << "Should not reach here inside lift_bdry. Exiting....." << "\n";
		exit(1);
	}
}

inline void compute_slopes(const std::vector<double> &var, std::vector<double> &slope, const int i)
{
	double deltaU_left, deltaU_right, slope1;

	deltaU_left =  var[i] - var[i-1];
	deltaU_right = var[i+1] - var[i];
		
	// Minmod limiter	
	slope1 = minmod(deltaU_left, deltaU_right);
	slope[i] = slope1;	
}

inline void state::compute_limited_slopes_consvar(int i)
{
	const auto &rho = dynamic.rho;
	const auto &rhou = dynamic.rhou;
	const auto &E = dynamic.E;
	
	// Compute slopes and store them 
	compute_slopes(rho , rho_slope, i);
	compute_slopes(rhou, rhou_slope, i);
	compute_slopes(E   , E_slope, i);	
}


inline void compute_left_and_right_states(const std::vector<double>& a_var, const int a_n, const int a_i, const std::vector<double>& a_slope,
                                      std::vector<double>& a_state_left, std::vector<double>& a_state_right)
{
    a_state_left[a_n]  = a_var[a_i-1] + 0.5*a_slope[a_i-1];
    a_state_right[a_n] = a_var[a_i] - 0.5*a_slope[a_i];
}

inline void state::compute_flux(int i)
{

	std::vector<double> QL(ND+2,0.0), QR(ND+2,0.0);

	const auto &rho = dynamic.rho;
	const auto &rhou = dynamic.rhou;
	const auto &E = dynamic.E;


	// Compute left and right states. ie. for all the conservative variables
	compute_left_and_right_states(rho,  0, i, rho_slope,  QL, QR);	
	compute_left_and_right_states(rhou, 1, i, rhou_slope, QL, QR);	
	compute_left_and_right_states(E   , 2, i, E_slope,    QL, QR);	

	std::array<double,3> flux = llf_flux(QL[0], QL[1], QL[2], QR[0], QR[1], QR[2]);

	flux_rho[i] = flux[0]; 
	flux_rhou[i] = flux[1]; 
	flux_E[i] = flux[2]; 
}

inline void state::write_solution(const int iter)
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
        		fprintf(output_file, "%0.15g %0.15g %0.15g %0.15g\n",x[i], dynamic.rho[i], dynamic.rhou[i], dynamic.E[i]);
    		}
			fflush(output_file);
		}	
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

    fclose(output_file);
}

inline void state::apply_bc(std::vector<double>& a_var)
{
	for(auto bdry: physical_boundaries)
	if(bdry == left){
		for(int i=0;i<=ng-1;i++){
			a_var[i] = a_var[ng];
		}
	}
	else if(bdry == right){	
		for(int i=sz+ng;i<=sz+2*ng-1;i++){
			a_var[i] = a_var[sz+ng-1];
		}
	}
}

inline void state::boundary_conditions(dynamic_state& a_dyn_state)
{
	apply_bc(a_dyn_state.rho);
	apply_bc(a_dyn_state.rhou);
	apply_bc(a_dyn_state.E);
}
}
