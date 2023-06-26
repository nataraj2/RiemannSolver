#include <iostream>
#include<vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <array>

//#include "Riemann.H"

int ND = 1;

double minmod(double a, double b);
template <class T1, class T2>
void ComputePrimitiveVariables(T1 &Q, T2 &prim);
void ComputeFluxes(const std::vector<double> &Q, std::array<double,3> &flux);
void ComputeQ_star(const std::vector<double> &Q, const double S_star, const double S,
                   const double rho, const double u, const double p,
                   std::vector<double> &Q_star);
void ComputeF_star(const std::array<double,3> &F, const std::vector<double> &Q,
                   const std::vector<double> &Q_star, const double S,
                   std::array<double,3> &F_star);
std::array<double,3> ComputeHLLCFlux(const double SL, const double SR,
                                     const std::array<double,3> &FL, const std::array<double,3> &FR,
                                     const double S_star, const std::array<double,3> &FL_star, const std::array<double,3> &FR_star);
void ComputeLeftAndRightStates(const std::vector<double>& a_var, const int a_comp, const int a_i, const std::vector<double>& a_slope,
                               std::vector<double>& a_dyn_state_left, std::vector<double>& a_dyn_state_right);
std::array<double,3> llf_flux(const double rho_l, const double rhou_l, const double E_l,
                                          const double rho_r, const double rhou_r, const double E_r);

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

template <class T1, class T2>
void ComputePrimitiveVariables(T1 &Q, T2 &prim)
{
	assert(Q[0] > 0.0);
	assert(Q[ND+1] > 0.0);
	prim[0] = Q[0];      // density
	prim[1] = Q[1]/Q[0]; // u-velocity
	prim[ND+1] = (Q[ND+1] - 0.5*prim[0]*prim[1]*prim[1])*(gamma_air-1.0); // pressure
}

void ComputeFluxes(const std::vector<double> &Q, std::array<double,3> &flux)
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

void ComputeF_star(const std::array<double,3> &F, const std::vector<double> &Q,
                   const std::vector<double> &Q_star, const double S,
                   std::array<double,3> &F_star)
{

	for(int n=0; n<=ND+1; n++){
		F_star[n] = F[n] + S*(Q_star[n] - Q[n]);
	}
}

std::array<double,3> ComputeHLLCFlux(const double SL, const double SR, 
					 const std::array<double,3> &FL, const std::array<double,3> &FR, 
					 const double S_star, const std::array<double,3> &FL_star, const std::array<double,3> &FR_star)
{
		std::array<double,3> flux;

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

std::array<double,3> llf_flux(const double rho_l, const double rhou_l, const double E_l,
    									  const double rho_r, const double rhou_r, const double E_r)
{

	double rhoL, uL, pL, cL, rhoR, uR, pR, cR, c_hat, u_hat;
	double SL, SR, S_star;
	std::vector<double> prim(ND+2,0.0);
	std::vector<double> QL_star(ND+2,0.0), QR_star(ND+2,0.0);

	std::array<double,3> FL, FR, FL_star, FR_star;

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

	std::array<double,3> flux = ComputeHLLCFlux(SL, SR, FL, FR, S_star, FL_star, FR_star);

	return flux;
}
