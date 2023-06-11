#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>
//#include <algorithm>

#include "Riemann.H"

int ND, nx, ny;
double xmin, xmax;

double rhoL, pL, uL;
double rhoR, pR, uR;
double gamma_air;
double dx;

std::vector<double> x;
std::vector<Cell> cell;
std::vector<Face> face;

void Initialize()
{
	
	ND = 1;
	nx = 501;
    xmin = -10.0, xmax = 10.0;

    rhoL = 1.0, pL = 100000.0, uL = 0.0;
    rhoR = 0.125, pR = 10000.0, uR = 0.0;
    //rhoR = 1.0, pR = 100000.0, uR = 0.0;
    gamma_air = 1.4;
	
	x.resize(nx+2);
	cell.resize(nx+2);
	face.resize(nx+3);
	
    dx = (xmax-xmin)/(nx-1);

    /*for(int i=0;i<=nx+1;i++){
        auto &icell = cell[i];
        x[i] = xmin -2.0*dx + i*dx;
        if(x[i] <= -0.01){
            icell.cons_var[0] = rhoL;
            icell.cons_var[1] = icell.cons_var[0]*uL;
            icell.cons_var[ND+1] = pL/(gamma_air-1.0) + 0.5*rhoL*uL*uL;

        }
        else if (x[i]>=0.01){
            icell.cons_var[0] = rhoR;
            icell.cons_var[1] = icell.cons_var[0]*uR;
            icell.cons_var[ND+1] = pR/(gamma_air-1.0) + 0.5*rhoR*uR*uR;
        }
		else{
			icell.cons_var[0] = rhoR*3.0;
			icell.cons_var[1] = icell.cons_var[0]*uR;
			icell.cons_var[ND+1] = 10.0*pR/(gamma_air-1.0) + 0.5*icell.cons_var[0]*uR*uR;
		}
	
    }*/

 	for(int i=0;i<=nx+1;i++){
        auto &icell = cell[i];
        x[i] = xmin -2.0*dx + i*dx;
        if(x[i] <= 0.0){
            icell.cons_var[0] = rhoL;
            icell.cons_var[1] = icell.cons_var[0]*uL;
            icell.cons_var[ND+1] = pL/(gamma_air-1.0) + 0.5*rhoL*uL*uL;

        }
        else if (x[i]>=0.0){
            icell.cons_var[0] = rhoR;
            icell.cons_var[1] = icell.cons_var[0]*uR;
            icell.cons_var[ND+1] = pR/(gamma_air-1.0) + 0.5*rhoR*uR*uR;
        }
    }


	

	for(int i=0;i<=nx+2;i++){
		face[i].flux.resize(ND+2);
	}

}

double minmod(double a, double b){

	if(a*b > 0 and std::fabs(a) < std::fabs(b)){
		return a;
	}
	else if(a*b > 0 and std::fabs(b) < std::fabs(a)){
		return b;
	}
	else if(a*b <= 0){
		return 0.0;
	}
	else{
		return 0.0;
	}
}	

double maxmod(double a, double b){

	if(a*b > 0 and std::fabs(a) > std::fabs(b)){
		return a;
	}
	else if(a*b > 0 and std::fabs(b) > std::fabs(a)){
		return b;
	}
	else if(a*b <= 0){
		return 0.0;
	}
	else{
		return 0.0;
	}
}	

void ComputeLimitedSlopes_Prim(std::vector<Cell> &cell)
{
	double deltaU_left, deltaU_right, deltaU_central, slope1, slope2;
	
	for(int i=0; i<=nx+1; i++){
		for(int n=0; n<=ND+1; n++){
			cell[i].slope[n] = 0.0;
		}
	}
	
	std::vector<double> prim(ND+1,0);
	
	for(int i=1; i<=nx; i++){
		for(int n=0; n<=ND+1; n++){
			deltaU_left =  cell[i].cons_var[n]   - cell[i-1].cons_var[n];
			deltaU_right = cell[i+1].cons_var[n] - cell[i].cons_var[n];
			/*slope1 = minmod(deltaU_left    , 2.0*deltaU_right);
			slope2 = minmod(2.0*deltaU_left, deltaU_right);
			cell[i].slope[n] = maxmod(slope1, slope2);*/

			// Minmod limiter	
			slope1 = minmod(deltaU_left    , deltaU_right);
			cell[i].slope[n] = slope1;	

			// Van Leer MC limiter
			//deltaU_central = cell[i+1].cons_var[n] - cell[i-1].cons_var[n];
			//slope1 = minmod(deltaU_central, 2.0*deltaU_left, 2.0*deltaU_right);
			//cell[i].slope[n] = slope1;	 
		}
	}	
}


void ComputeLimitedSlopes_ConsVar(std::vector<Cell> &cell)
{
	double deltaU_left, deltaU_right, deltaU_central, slope1, slope2;
	
	for(int i=0; i<=nx+1; i++){
		for(int n=0; n<=ND+1; n++){
			cell[i].slope[n] = 0.0;
		}
	}
		
	for(int i=1; i<=nx; i++){
		for(int n=0; n<=ND+1; n++){
			deltaU_left =  cell[i].cons_var[n]   - cell[i-1].cons_var[n];
			deltaU_right = cell[i+1].cons_var[n] - cell[i].cons_var[n];
			/*slope1 = minmod(deltaU_left    , 2.0*deltaU_right);
			slope2 = minmod(2.0*deltaU_left, deltaU_right);
			cell[i].slope[n] = maxmod(slope1, slope2);*/

			// Minmod limiter	
			slope1 = minmod(deltaU_left    , deltaU_right);
			cell[i].slope[n] = slope1;	

			// Van Leer MC limiter
			//deltaU_central = cell[i+1].cons_var[n] - cell[i-1].cons_var[n];
			//slope1 = minmod(deltaU_central, 2.0*deltaU_left, 2.0*deltaU_right);
			//cell[i].slope[n] = slope1;	 
		}
	}	
}

template <class T>
void ComputePrimVariables(T &Q, std::vector<double> &prim)
{
	assert(Q[0] > 0.0);
	assert(Q[ND+1] > 0.0);
	prim[0] = Q[1]/Q[0]; // u-velocity
	prim[1] = (Q[ND+1] - 0.5*Q[0]*prim[0]*prim[0])*(gamma_air-1.0); // pressure
}

void ComputeFluxes(const std::vector<double> &Q, std::vector<double> &flux)
{
	
	double uvel = Q[1]/Q[0];
	double pressure = (Q[ND+1] - 0.5*Q[0]*uvel*uvel)*(gamma_air-1.0);
	assert(pressure > 0.0);
	flux[0] = Q[1];
	flux[1] = pressure + Q[0]*uvel*uvel;
	flux[2] = (Q[ND+1] + pressure)*uvel;

}

void ComputeQ_star(const std::vector<double> &Q, const double S_star, const double S, 
				   const double rho, const double u, const double p, 
				   std::vector<double> &Q_star)
{

	double fac = rho*(S-u)/(S-S_star); 

	Q_star[0] = fac;
	Q_star[1] = fac*S_star;
	Q_star[ND+1] = fac*(Q[ND+1]/Q[0] + (S_star - u)*(S_star + p/(rho*(S-u))));	

}

void ComputeF_star(const std::vector<double> &F, const std::vector<double> &Q,
                   const std::vector<double> &Q_star, const double S,
                   std::vector<double> &F_star)
{

	for(int n=0; n<=ND+1; n++){
		F_star[n] = F[n] + S*(Q_star[n] - Q[n]);
	}
}

void ComputeHLLCFlux(const double SL, const double SR, 
					 const std::vector<double> &FL, const std::vector<double> &FR, 
					 const double S_star, const std::vector<double> &FL_star, const std::vector<double> &FR_star, 
					 Face &iface)
{
		if(0.0 <= SL){
			iface.flux = FL;
		}
		else if(SL<= 0.0 and 0.0 <=S_star){
			iface.flux = FL_star;
		}
		else if(S_star <=0.0 and 0.0 <= SR){
			iface.flux = FR_star;
		}
		else if(0.0 >= SR){
			iface.flux = FR;
		}
		else{
			std::cout << "Reached HLLC flux unknown" << "\n";
			exit(1);
		}
}

void ComputeHLLCFluxOnFace(const std::vector<double> &QL, const std::vector<double> &QR, Face &iface)
{

	double rhoL, uL, pL, cL, rhoR, uR, pR, cR, c_hat, u_hat;
	double SL, SR, S_star;
	std::vector<double> prim(ND+1,0.0), FL(ND+2,0.0), FR(ND+2,0.0);
	std::vector<double> FL_star(ND+2,0.0), FR_star(ND+2,0.0), QL_star(ND+2,0.0), QR_star(ND+2,0.0);
	
	ComputePrimVariables(QL,prim);

	rhoL = QL[0];
	uL   = prim[0];
	pL   = prim[1];
	
	ComputePrimVariables(QR,prim);

	rhoR = QR[0];
	uR   = prim[0];
	pR   = prim[1];

	cL = sqrt(gamma_air*pL/rhoL);
	cR = sqrt(gamma_air*pR/rhoR);


	u_hat = (uL*sqrt(rhoL) + uR*sqrt(rhoR))/(sqrt(rhoL) + sqrt(rhoR));
	c_hat = sqrt((std::pow(cL,2)*std::sqrt(rhoL)+std::pow(cR,2)*std::sqrt(rhoR))/(std::sqrt(rhoL)+std::sqrt(rhoR)) +
				 std::sqrt(rhoL*rhoR)/(2.0*pow((std::sqrt(rhoL)+std::sqrt(rhoR)),2))*std::pow(uR-uL,2));

	SL = std::min(uL-cL,u_hat-c_hat); 
	SR = std::min(uR+cR,u_hat+c_hat); 

	ComputeFluxes(QL, FL);
	ComputeFluxes(QR, FR);

	S_star = (pR- pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR))/(rhoL*(SL-uL) - rhoR*(SR-uR));

	ComputeQ_star(QL, S_star, SL, rhoL, uL, pL, QL_star);
	ComputeQ_star(QR, S_star, SR, rhoR, uR, pR, QR_star);

	ComputeF_star(FL, QL, QL_star, SL, FL_star);		
	ComputeF_star(FR, QR, QR_star, SR, FR_star);		
			
	ComputeHLLCFlux(SL, SR, FL, FR, S_star, FL_star, FR_star, iface);
}

void ComputeLeftAndRightStatesAndComputeFlux(std::vector<Cell> &cell, std::vector<Face> &face)
{

	std::vector<double> QL(ND+2,0.0), QR(ND+2,0.0);

	// Face loop over all faces
	for(int i=2;i<=nx; i++){
		// Compute left and right states. ie. all the conservative variables
		for(int n=0;n<=ND+1;n++){
			QL[n] = cell[i-1].cons_var[n] + 0.5*cell[i-1].slope[n];  
			QR[n] = cell[i].cons_var[n] - 0.5*cell[i].slope[n];
		}
		//printf("%d, %g, %g, %g\n", i, QL[0], QL[1], QL[2]);
		//printf("%d, %g, %g, %g\n", i, QR[0], QR[1], QR[2]);
		ComputeHLLCFluxOnFace(QL, QR, face[i]);	
	}
}

void BoundaryConditions(std::vector<Cell> &cell)
{

	std::copy(std::begin(cell[2].cons_var), std::end(cell[2].cons_var), std::begin(cell[0].cons_var));
	std::copy(std::begin(cell[2].cons_var), std::end(cell[2].cons_var), std::begin(cell[1].cons_var));
	std::copy(std::begin(cell[nx-1].cons_var), std::end(cell[nx-1].cons_var), std::begin(cell[nx].cons_var));
	std::copy(std::begin(cell[nx-1].cons_var), std::end(cell[nx-1].cons_var), std::begin(cell[nx+1].cons_var));

}

double ComputeTimeStep(const std::vector<Cell> &cell, const double CFL)
{
	std::vector<double> prim(ND+1,0.0);
	std::vector<double> uvel(nx+2,0.0);
	std::vector<double> spdsnd(nx+2,0.0);
	
	for(int i=0;i<=nx+1;i++){
		ComputePrimVariables(cell[i].cons_var, prim);
		uvel[i]    = prim[0];
		spdsnd[i] = sqrt(gamma_air*prim[1]/cell[i].cons_var[0]);
	}

	auto uvel_max   = std::max_element(uvel.begin(), uvel.end());
	auto spdsnd_max = std::max_element(spdsnd.begin(), spdsnd.end());
	
	double dt = CFL*dx/std::max(*uvel_max, *spdsnd_max);
	
	return dt;
}

void WriteSolution(const std::vector<Cell> &cell, const int iter)
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
        auto &icell = cell[i];
        fprintf(output_file, "%g %g %g %g\n",x[i], icell.cons_var[0], icell.cons_var[1], icell.cons_var[ND+1]);
    }

    fclose(output_file);

}

