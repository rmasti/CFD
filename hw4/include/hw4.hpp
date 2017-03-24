/*
 * Comp. Fluid Dynamcis
 * Written By: Robert Masti
 * 1/27/2017
 * This is a header file for nozzlem.cpp which will contain calls from the nozzlef.cpp file which will contain universal constants as well as function prototypes, and even structure declarations.
 */
#ifndef hw4_H_
#define hw4_H_
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "math.h"
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;
using namespace Eigen;

#define N 256
#define xmax 1.5
#define xmin -1.5
#define xmax_dom 1.0
#define xmin_dom -1.0
#define num_ghost_cells 3
#define R 287.058
#define kappa2 (1.0/3.0)
#define kappa4 (1.0/48.0)
#define dx (xmax-xmin)/(N)
#define localdt false
#define nmax 10000000
#define mymax(a,b) ((a>b)?a:b)
#define mymin(a,b) ((a<b)?a:b)

#define rhoid  0
#define uid 1
#define vid 2
#define pid 3
#define rhouid 1
#define rhovid 2
#define rhoetid 3
#define frhouid 0
#define frhouuid 1
#define frhovid 2
#define frhouhtid 3
#define neq 4

// CREATE STRUCTURE DEFINITION
struct constants
{
  double p0;
  double T0;
  double A_t;
  double tol;
  double gamma;
  bool outflow;
  double pb; //back pressure
  double cfl;
};

////////////// Function Prototypes ////////////////////
// Main
void output_file_headers();

void quasi1Dnozzle(constants C);

// Function File
double A_x(double x);

double M_xinitial(double x);

double dAdx(double x);

void set_geometry(MatrixXd& xcenter, MatrixXd& xinterface);   

void initialize(MatrixXd* V, MatrixXd &Mc, constants &C);

double compute_soundspeed( double gamma,  double p, double rho);

void Mtoprim(double& V1, double& V2, double& V3, double& V4, 
    double& M, constants& C);

double compute_soundspeed(double gamma, double p, double rho);

void extrapolate_to_ghost(MatrixXd* Varr);

void set_boundary_conditions( MatrixXd* V, constants C);

void inflow_boundary_condition(MatrixXd* V, constants C);

void outflow_boundary_condition(MatrixXd* V, constants C);

void primtocons(MatrixXd* U, MatrixXd* V, constants C);     

void compute_lambda(MatrixXd& Lamda_mcenter, MatrixXd* V, constants C);

void reconstruct(MatrixXd* Uinterface, MatrixXd& Lambda_minterface, 
    MatrixXd* U, MatrixXd& Lambda_mcenter); 

void compute_F_flux(MatrixXd* F, MatrixXd* U, constants C);

void artificial_viscosity(MatrixXd* d, MatrixXd* V, 
    MatrixXd* U, MatrixXd& Lambda_minterface);

void compute_nu(MatrixXd& nu, MatrixXd* V);

void compute_source(MatrixXd* S, MatrixXd* V, MatrixXd& xc);

void compute_residual(MatrixXd* Res, MatrixXd* S, MatrixXd* F,
    MatrixXd* d, MatrixXd& Ai);

void iteration(MatrixXd* U, double& timestep, MatrixXd* Res, 
    MatrixXd& Ac, MatrixXd& Lambda_mcenter, constants C);

void constoprim(MatrixXd* V, MatrixXd* U, constants C);

void exactsol(MatrixXd* V, MatrixXd& Ac, constants C);

void isentropic_exact(constants C);

void compute_norms(MatrixXd& L2norm, MatrixXd* Res);

void write_solution(FILE* &file, MatrixXd& xc, MatrixXd& Ac,
    MatrixXd* V, MatrixXd* U, constants C);

#endif
