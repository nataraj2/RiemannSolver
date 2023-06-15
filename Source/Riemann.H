class Cell{
public:
    double cons_var[3], prim_var[3];
	double slope[3];

	//Cell(double rho_, double rhoU_, double rhoE_): 
	//	rho(rho_), rhoU(rhoU_), rhoE(rhoE_) {};
};

class Face{
	public:
		std::vector<double> flux;
};

extern int ND;
extern int nx, ny;
extern double xmin, xmax;
extern double gamma_air;
extern double dx, dt, t_end;

extern std::vector<double> x;
extern std::vector<Cell> cell;
extern std::vector<Face> face;

void Initialize();
void ComputeLeftAndRightCells(const int i, Cell &cell_left, Cell &cell_right);
double minmod(double a, double b);
double maxmod(double a, double b);
void ComputeLimitedSlopes_ConsVar(std::vector<Cell> &cell);
void ComputeLimitedSlopes_PrimVar(std::vector<Cell> &cell);
template <class T1, class T2>
void ComputePrimitiveVariables(T1 &Q, T2 &prim);
void ComputeFluxes(const std::vector<double> &Q, std::vector<double> &flux);
void ComputeQ_star(const std::vector<double> &Q, const double S_star, const double S,
                   const double rho, const double u, const double p,
                   std::vector<double> &Q_star);
void ComputeF_star(const std::vector<double> &F, const std::vector<double> &Q, 
				   const std::vector<double> &Q_star, const double S, 
				   std::vector<double> &F_star);
std::vector<double> ComputeHLLCFlux(const double SL, const double SR,
                     const std::vector<double> &FL, const std::vector<double> &FR,
                     const double S_star, const std::vector<double> &FL_star, const std::vector<double> &FR_star);
void ComputeLeftAndRightConsVarStatesAndComputeFlux(std::vector<Cell> &cell, std::vector<Face> &face);
void ComputeLeftAndRightPrimVarStatesAndComputeFlux(std::vector<Cell> &cell, std::vector<Face> &face);
std::vector<double> ComputeHLLCFluxOnFace(const double rho_l, const double rhou_l, const double E_l,
   		 								  const double rho_r, const double rhou_r, const double E_r);
void BoundaryConditions(std::vector<Cell> &cell);
double ComputeTimeStep(const std::vector<Cell> &cell, const double CFL);
void WriteSolution(const std::vector<Cell> &cell, const int iter);
