/*
 * Comp. Fluid Dynamcis
 * Written By: Robert Masti
 * 1/27/2017
 * This is a header file for nozzlem.cpp which will contain calls from the nozzlef.cpp file which will contain universal constants as well as function prototypes, and even structure declarations.
 */
#ifndef nozzle_H_
#define nozzle_H_
#define N 150
#define xmax 1.5
#define xmin -1.5
#define xmax_dom 1.0
#define xmin_dom -1.0
#define num_ghost_cells 3
#define R 287.0
#define kappa4 (1.0/64.0)
#define kappa2 (1.0/4.0)
#define cfl 0.1
#define dx (xmax-xmin)/(N)
#define localdt false
#define nmax 50000
/////////////////////////////////////////////////////////////////////////
///////////////////////// STRUCTURE DEFINITIONS /////////////////////////
/////////////////////////////////////////////////////////////////////////

// prim_exact is a structure that containing structure types double 
struct prim_exact
{
  double rho;
  double u;
  double p;
  double M;
};

struct primvar
{
  double rho;
  double u;
  double v;
  double p;
};

// constants contains doubles and booleans
struct constants
{
  double p0;
  double T0;
  double A_t;
  bool cond; //0 for subsonic vs //1 for 
  double tol;
  double gamma;
  bool outflow;
  double pb; //back pressure
};
// create consvar structure uses type double
struct consvar
{
  double rho;
  double rhou;
  double rhov;
  double rhoet;
};
struct fluxes
{
  double rhou;
  double rhouu_and_p;
  double rhouht;
};

/////////////////////////////////////////////////////////////////////////
///////////////////////// FUNCTIONS PROTOTYPES //////////////////////////
/////////////////////////////////////////////////////////////////////////
primvar constoprim(consvar U, constants C);
consvar primtocons(primvar V, constants C);
prim_exact exactsol(double A_x, constants C);
primvar Mtoprim(double M, constants C);

double A_x(double xcoord);
double dAdx(double x);
double M_x(double x);

void isentropic(constants C);

void isentropic_exact(constants C);

void initialize(std::vector<double> &M, std::vector<consvar> &U, constants C);

void set_geometry(std::vector<double> &Aarr, std::vector<double> & xinterface, std::vector<double> & xcenter, std::vector<double> &Marr);

void output_file_headers();

void write_out(FILE* &fp2, std::vector<double> const &Aarr, std::vector<double> const &XCarr, std::vector<consvar> const &U, constants C);

void set_boundary_cond(std::vector<consvar> &U, constants C);

void compute_fluxes(std::vector<fluxes> &F, std::vector<consvar> const &U, constants C);

void reconstruct_U(std::vector<consvar> &U_avg, std::vector<consvar> const &U);

fluxes fluxcalc( consvar const &U, constants C);

void iteration_step(std::vector<fluxes> &F, std::vector<consvar> &Uold, std::vector<consvar> &Unew, std::vector<primvar> &Vold, std::vector<primvar> &Vnew, std::vector<double> const &XCarr, std::vector<double> const &Xarr, std::vector<double> &Marr, std::vector<consvar> &Resarr, constants C);

void iteration_step(std::vector<fluxes> &F, std::vector<consvar> &Uold, std::vector<consvar> &Unew, std::vector<double> const &XCarr, std::vector<double> const &Xarr, std::vector<consvar> &Resarr, constants C);

double compute_timestep(std::vector<consvar> const &Uold, int i, constants C);

double compute_volume(std::vector<double> const &Xarr, int i, std::vector<double> &ALR);

double primtoM(primvar V, constants C);

void compute_residuals(std::vector<consvar> &Resarr, std::vector<double> &Res, std::vector<double> &Linfnorm, std::vector<double> &L1norm, std::vector<double> &L2norm, std::vector<consvar> const &Uold, std::vector<consvar> const &Unew);

void artificial_viscosity(std::vector<fluxes> &dvec, std::vector<consvar> const &U, constants C);

fluxes compute_dflux( double const &epsilon2, double const &epsilon4, double const &lambda, int const   &i, std::vector<consvar> const &U);

void extrapolate_to_ghost(std::vector<consvar> &Uarr, constants C);

#endif

