/*
 * Comp. Fluid Dynamcis
 * Written By: Robert Masti
 * 1/27/2017
 * This is a header file for nozzlem.cpp which will contain calls from the nozzlef.cpp file which will contain universal constants as well as function prototypes, and even structure declarations.
 */
#ifndef nozzle_H_
#define nozzle_H_
#define N 300
#define xmax 1.5
#define xmin -1.5
#define xmax_dom 1.0
#define xmin_dom -1.0


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
  double dx;
  double nmax;
};
// create consvar structure uses type double
struct consvar
{
  double rho;
  double rhou;
  double rhov;
  double rhoet;
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

void isentropicExact(constants C);
void initialize(std::vector<primvar> &V, std::vector<double> const &M, constants C);
void set_geometry(std::vector<double> &Aarr, std::vector<double> &Xarr, 
    std::vector<double> &dAdxarr, std::vector<double> &Marr);

#endif
