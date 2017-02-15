/*
 * Computational Fluid Dynamics
 * Written By: Robert Masti
 * 1/27/2017
 * This file contains all function declarations that will be called in hw1_m.cpp which is the 
 * main cpp file. This is the function file, which has its function prototypes stored in the 
 * header file.
 */
#include <iostream>
#include <cmath>
#include "nozzle.h" //structure templates and func prototypes
#include <fstream>
#include <iomanip>
using namespace std;
// This function converts a conservative variable vector into a primitive variable vector
primvar constoprim(consvar U, constants C)
{
  primvar out;
  out.rho = U.rho;
  double invrho = 1.0/U.rho;
  out.u = U.rhou*invrho;
  out.v = U.rhov*invrho;
  out.p = (C.gamma-1.0)*(U.rhoet - 0.5*invrho*U.rhou*U.rhou);
  return out;
}

// Compute the area and derivative of area it takes in double x and return double A
double A_x(double x)
{
  using namespace std;
  if (x < xmin || x > xmax) cout << "ERROR: x value is greater then specified domain. Change xmax in header" << endl;
  
  if (x < xmin_dom || x > xmax_dom) return 1.0;
  else return 0.2 + 0.4*(1.0 + sin(M_PI*(x-0.5)));//function that governs Area for nozzle
}
double dAdx(double x)
{
  if (x < -1.0 || x > 1.0) return 0.0;
  else return 0.4*M_PI*cos(M_PI*(x-0.5));//derivative of func above
}

//compute exact solution
prim_exact exactsol(double A_x, constants stagpTAcond)//constants is a datastruct defined in the header
{
  /*************************************************************************************/
  /********************************* Local Definitions *********************************/
  /*************************************************************************************/
  using namespace std;
  prim_exact answer;
  double A_t = stagpTAcond.A_t;//m^2
  double A_bar = A_x/A_t;//unitless
  double p0 = (stagpTAcond.p0)*1000.0;//kPa->Pa
  double T0 = stagpTAcond.T0;//K
  double gam = 1.4;//unitless
  double tol = stagpTAcond.tol;//unitless
  double phi;//unitless
  double F = 2.0;//unitless
  double dFdM;//unitless
  double Mnew;//unitless
  double dM=100000.0;//unitless
  //must compute the mach number and then with that one can compute the state variables
  double M;

  /*************************************************************************************/
  /**************************** Create the initial guess of Mach ***********************/
  /*************************************************************************************/
  if (stagpTAcond.cond == true) M=0.1;
  else M=5.0;

  /*************************************************************************************/
  /******************************* Newton Raphson Loop *********************************/
  /*************************************************************************************/
  while (abs(dM) > tol)
  {
    phi = (2.0/(gam + 1.0)) * (1.0 + ((gam - 1.0)/2.0) * pow(M,2.0));//unitless
    F = pow(phi,(gam + 1.0)/(gam - 1.0)) - pow(A_bar,2.0) * pow(M,2.0);//unitless
    dFdM = 2.0 * M * (pow(phi,2.0/(gam - 1.0)) - pow(A_bar,2.0));//unitless
    Mnew = M - (F/dFdM); //unitless 
    dM = (Mnew-M);//Used to check convergence
    M = Mnew;
  }
  //return the primary variables V once the Mach number is known.
  /*************************************************************************************/
  /***************************** Compute State Variables *******************************/
  /*************************************************************************************/
  double psi = 1.0 + ((gam - 1.0)/2.0) * pow(M,2.0);//unitless
  double T = T0/psi;// K
  double R = 287.0;//J/(kg K)
  answer.p = (p0/pow(psi,gam/(gam - 1.0)));//Pa
  answer.rho = answer.p/(R * T);
  answer.u = M * sqrt(gam * R * T);
  answer.M = M;
  return answer;//returns a prim_exact data struct defined in the header file.
}
