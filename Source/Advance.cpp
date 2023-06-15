#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "Riemann.H"

int ND, nx, ny;
double xmin, xmax;
double gamma_air;
double dx, dt, t_end;

std::vector<double> x;
std::vector<Cell> cell;
std::vector<Face> face;


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

void ComputeLimitedSlopes_PrimVar(std::vector<Cell> &cell)
{
	double deltaU_left, deltaU_right, deltaU_central, slope1, slope2;
	
	for(int i=0; i<=nx+1; i++){
		for(int n=0; n<=ND+1; n++){
			cell[i].slope[n] = 0.0;
		}
	}
	
	for(int i=0; i<=nx+1; i++){
		ComputePrimitiveVariables(cell[i].cons_var, cell[i].prim_var);
	}


	for(int i=0; i<=nx+1; i++){
		for(int n=0; n<=ND+1; n++){
			deltaU_left =  cell[i].prim_var[n]   - cell[i-1].prim_var[n];
			deltaU_right = cell[i+1].prim_var[n] - cell[i].prim_var[n];

			// Superbee limiter
			//slope1 = minmod(deltaU_left    , 2.0*deltaU_right);
			//slope2 = minmod(2.0*deltaU_left, deltaU_right);
			//cell[i].slope[n] = maxmod(slope1, slope2);

			// Minmod limiter	
			slope1 = minmod(deltaU_left, deltaU_right);
			cell[i].slope[n] = slope1;	
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

template <class T1, class T2>
void ComputePrimitiveVariables(T1 &Q, T2 &prim)
{
	assert(Q[0] > 0.0);
	assert(Q[ND+1] > 0.0);
	prim[0] = Q[0];      // density
	prim[1] = Q[1]/Q[0]; // u-velocity
	prim[ND+1] = (Q[ND+1] - 0.5*prim[0]*prim[1]*prim[1])*(gamma_air-1.0); // pressure
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

std::vector<double> ComputeHLLCFlux(const double SL, const double SR, 
					 const std::vector<double> &FL, const std::vector<double> &FR, 
					 const double S_star, const std::vector<double> &FL_star, const std::vector<double> &FR_star)
{
		std::vector<double> flux(ND+2,0.0);

		if(0.0 <= SL){
			flux = FL;
		}
		else if(SL<= 0.0 and 0.0 <=S_star){
			flux = FL_star;
		}
		else if(S_star <=0.0 and 0.0 <= SR){
			flux = FR_star;
		}
		else if(0.0 >= SR){
			flux = FR;
		}
		else{
			std::cout << "Reached HLLC flux unknown" << "\n";
			exit(1);
		}

		return flux;
}

std::vector<double> ComputeHLLCFluxOnFace(const double rho_l, const double rhou_l, const double E_l,
    const double rho_r, const double rhou_r, const double E_r)
{

	double rhoL, uL, pL, cL, rhoR, uR, pR, cR, c_hat, u_hat;
	double SL, SR, S_star;
	std::vector<double> prim(ND+2,0.0), FL(ND+2,0.0), FR(ND+2,0.0);
	std::vector<double> FL_star(ND+2,0.0), FR_star(ND+2,0.0), QL_star(ND+2,0.0), QR_star(ND+2,0.0);

	std::vector<double> QL(ND+2,0.0), QR(ND+2,0.0);

	QL[0] = rho_l, QL[1] = rhou_l, QL[2] = E_l;
	QR[0] = rho_r, QR[1] = rhou_r, QR[2] = E_r;
	
	ComputePrimitiveVariables(QL,prim);

	rhoL = prim[0];
	uL   = prim[1];
	pL   = prim[2];
	
	ComputePrimitiveVariables(QR,prim);

	rhoR = prim[0];
	uR   = prim[1];
	pR   = prim[2];

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

	std::vector<double> flux = ComputeHLLCFlux(SL, SR, FL, FR, S_star, FL_star, FR_star);

	return flux;
}

void ComputeLeftAndRightConsVarStatesAndComputeFlux(std::vector<Cell> &cell, std::vector<Face> &face)
{

	std::vector<double> QL(ND+2,0.0), QR(ND+2,0.0);

	// Face loop over all faces
	for(int i=2;i<=nx; i++){
		// Compute left and right states. ie. all the conservative variables
		for(int n=0;n<=ND+1;n++){
			QL[n] = cell[i-1].cons_var[n] + 0.5*cell[i-1].slope[n];  
			QR[n] = cell[i].cons_var[n] - 0.5*cell[i].slope[n];
		}

		
		//ComputeHLLCFluxOnFace(QL[0], QL[1], QL[2], QR[0], QR[1], QR[2], face[i]);	
	}
}

void ComputeConservativeVariables(const std::vector<double> &P, std::vector<double> &Q) 
{
	
	Q[0] = P[0];
	Q[1] = P[0]*P[1];
	Q[ND+1] = P[ND+1]/(gamma_air-1) + 0.5*P[0]*P[1]*P[1];
	
}

void ComputeLeftAndRightPrimVarStatesAndComputeFlux(std::vector<Cell> &cell, std::vector<Face> &face)
{

	std::vector<double> PL(ND+2,0.0), PR(ND+2,0.0);
	std::vector<double> QL(ND+2,0.0), QR(ND+2,0.0);

	// Face loop over all faces
	for(int i=2;i<=nx; i++){
		// Compute left and right states. ie. all the conservative variables
		for(int n=0;n<=ND+1;n++){
			PL[n] = cell[i-1].prim_var[n] + 0.5*cell[i-1].slope[n];  
			PR[n] = cell[i].prim_var[n] - 0.5*cell[i].slope[n];
		}
		ComputeConservativeVariables(PL, QL);
		ComputeConservativeVariables(PR, QR);

		face[i].flux = ComputeHLLCFluxOnFace(QL[0], QL[1], QL[2], QR[0], QR[1], QR[2]);	
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
	std::vector<double> prim(ND+2,0.0);
	std::vector<double> uvel(nx+2,0.0);
	std::vector<double> spdsnd(nx+2,0.0);
	
	for(int i=0;i<=nx+1;i++){
		ComputePrimitiveVariables(cell[i].cons_var, prim);
		uvel[i]    = prim[1];
		spdsnd[i] = sqrt(gamma_air*prim[2]/prim[0]);
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

