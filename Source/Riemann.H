class Cell{
public:
    double cons_var[3], prim_var[3];
	double slope_cons_var[3], slope_prim_var[3];

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
extern double dx, dt;

extern std::vector<double> x;
extern std::vector<Cell> cell;
extern std::vector<Face> face;

void Initialize();
void ComputeLeftAndRightCells(const int i, Cell &cell_left, Cell &cell_right);
double minmod(double a, double b);
double maxmod(double a, double b);
void ComputeLimitedSlopes_ConsVar(std::vector<Cell> &cell);
template <class T>
void ComputePrimitiveVariables(T &Q, std::vector<double> &prim);
void ComputeFluxes(const std::vector<double> &Q, std::vector<double> &flux);
void ComputeQ_star(const std::vector<double> &Q, const double S_star, const double S,
                   const double rho, const double u, const double p,
                   std::vector<double> &Q_star);
void ComputeF_star(const std::vector<double> &F, const std::vector<double> &Q, 
				   const std::vector<double> &Q_star, const double S, 
				   std::vector<double> &F_star);
void ComputeHLLCFlux(const double SL, const double SR,
                     const std::vector<double> &FL, const std::vector<double> &FR,
                     const double S_star, const std::vector<double> &FL_star, const std::vector<double> &FR_star,
                     Face &iface);
void ComputeLeftAndRightStatesAndComputeFlux(std::vector<Cell> &cell, std::vector<Face> &face);
void ComputeHLLCFluxOnFace(const std::vector<double> &QL, const std::vector<double> &QR, Face &iface);
void BoundaryConditions(std::vector<Cell> &cell);
double ComputeTimeStep(const std::vector<Cell> &cell, const double CFL);
void WriteSolution(const std::vector<Cell> &cell, const int iter);
