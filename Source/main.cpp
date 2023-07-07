#include <iostream>
#include <string>
#include <vector>

#include "globals.hpp"
#include "ics.hpp"
#include "state.hpp"
#include "mpi.h"

using namespace Euler;

int rank_me, rank_n;

struct BoundaryWrapper {
    BoundaryWrapper()=default;
    void initialize(uint64_t actor_id, MPI_Request* send_req, MPI_Request* recv_req) {

        id = actor_id;

        MPI_Send_init(send_buffer.data(),
                      ng*3,
                      MPI_DOUBLE,
                      id,
                      0,
                      MPI_COMM_WORLD,
                      send_req
            );

        MPI_Recv_init(receive_buffer.data(),
                      ng*3,
                      MPI_DOUBLE,
                      id,
                      0,
                      MPI_COMM_WORLD,
                      recv_req
            );

    }

    uint64_t id;

	std::array<std::array<double,3>,ng> send_buffer, receive_buffer;
};


int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_me);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_n);

	double start_time, end_time, elapsed_time;

    start_time = MPI_Wtime();

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

	bool has_left_bdry, has_right_bdry;
	std::vector<BoundaryWrapper> distributed_boundaries;
	std::vector<MPI_Request> send_request, receive_request;	

	if(rank_me==0){
		has_left_bdry = false;
		has_right_bdry = true;
	}
	else if(rank_me==rank_n-1){
		has_left_bdry = true;
        has_right_bdry = false;
	}
	else{
		has_left_bdry = true;
        has_right_bdry = true;
	}


	size_t n_bdries = (bool)has_left_bdry + (bool)has_right_bdry;
    send_request.resize(n_bdries);
    receive_request.resize(n_bdries);
    distributed_boundaries.resize(n_bdries);
    if ( has_left_bdry ) {
        distributed_boundaries[0].initialize(rank_me-1,
                                             &send_request[0],
                                             &receive_request[0]);
		distributed_boundaries[0].id = rank_me-1;
    }
    if ( has_right_bdry ) {
        distributed_boundaries[(bool)has_left_bdry].initialize(rank_me+1,
                                                              &send_request[(bool)has_left_bdry],
                                                              &receive_request[(bool)has_left_bdry]);
		distributed_boundaries[(bool)has_left_bdry].id = rank_me+1;
    }

	int id = rank_me;

	state_.allocate_arrays();
	
	// Time advance
	for(int iter=1; iter<=niter; iter++){
		for ( auto& dist_bdry : distributed_boundaries ){
        	auto msg_outgoing = state_.make_message(id, dist_bdry.id);
        	msg_outgoing.write_to_array(dist_bdry.send_buffer);
    	}

    	MPI_Startall(send_request.size(), send_request.data());
    	MPI_Startall(receive_request.size(), receive_request.data());
			
		state_.build_u_next(state_.dynamic, rank_me, rank_n);

		MPI_Waitall(send_request.size(), send_request.data(), MPI_STATUSES_IGNORE);
    	MPI_Waitall(receive_request.size(), receive_request.data(), MPI_STATUSES_IGNORE);

		for ( auto& dist_bdry : distributed_boundaries ) {
			state_.lift_bdry(id, dist_bdry.id, state_.dynamic, dt*state_.process_message(id, dist_bdry.id,
                                                msg_t(dist_bdry.receive_buffer)));
		}

		if(iter%freq==0){
			state_.write_solution(iter/freq);
		}

		time = time + dt;
		if(rank_me == 0){
			std::cout << "Iteration is " << iter << " dt = " << dt << " " << "time = " << time << "\n";
		}
	}
	end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;

    printf("Process %d: Elapsed time = %.6f seconds\n", rank_me, elapsed_time);

	MPI_Barrier(MPI_COMM_WORLD);
	return 0;

}

