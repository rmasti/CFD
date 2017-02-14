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

/////////////////////////////////////////////////////////////////////////
///////////////////////// STRUCTURE DEFINITIONS /////////////////////////
/////////////////////////////////////////////////////////////////////////

// primvar is a structure that containing structure types double 
struct primvar
{
  double rho;
  double u;
  double p;
  double M;
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
  double rhoet;
};

/////////////////////////////////////////////////////////////////////////
///////////////////////// FUNCTIONS PROTOTYPES //////////////////////////
/////////////////////////////////////////////////////////////////////////
primvar constoprim(consvar U, constants C);
primvar exactsol(double A_x, constants stagpTAcond);
double A_x(double xcoord);
double dAdx(double x);
void isentropicExact(constants consts);

#endif
