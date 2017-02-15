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
#include <fstream>
#include <iomanip>
#include <vector>

#include "nozzle.h" //structure templates and func prototypes

using namespace std;

void initialize(vector<primvar> &V)
{
  for (int i=0 ; i<N; i++)
  {
    primvar temp;
    temp.rho = 1.0*i;
    temp.u = 2.0*i;
    temp.v = 0.0;
    temp.p = 3.0*i;
    V.push_back(temp);
  }
}

/*******************************************************************************/
/********************************* SET GEOMETRY ********************************/
/*******************************************************************************/
void set_geometry(vector<double> &Aarr, vector<double> &Xarr, 
    vector<double> &dAdxarr, vector<double> &Marr)
{
  double dx = (xmax-xmin)/N;
  for (int i=0; i<N ; i++)
  {
    double x = xmin + dx*i;
    Xarr.push_back(x);
    Aarr.push_back(A_x(x));
    dAdxarr.push_back(dAdx(x));
    Marr.push_back(M_x(x));
  }
}
// Compute the area and derivative of area it takes in double x and return double A. It also computes initial mach number guess.
double A_x(double x)
{
  if (x < xmin_dom || x > xmax_dom) return 1.0;
  else return 0.2 + 0.4*(1.0 + sin(M_PI*(x-0.5)));//function that governs Area for nozzle
}
double dAdx(double x)
{
  if (x < -1.0 || x > 1.0) return 0.0;
  else return 0.4*M_PI*cos(M_PI*(x-0.5));//derivative of func above
}
double M_x(double x)
{
  double slope = (2.0 - 0.1)/(xmax_dom - xmin_dom);
  double b = 0.1-slope*xmin_dom;
  if (x < xmin_dom) return slope*xmin_dom+b;
  if (x > xmax_dom ) return slope*xmax_dom+b;
  else return slope*x+b;
}
/*******************************************************************************/
/****************************** END SET GEOMETRY *******************************/
/*******************************************************************************/


/*******************************************************************************/
/********************************* CONVERSIONS *********************************/
/*******************************************************************************/

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
// This function does vice versa
consvar primtocons(primvar V, constants C)
{
  consvar out;
  out.rho = V.rho;
  out.rhou = V.rho*V.u;
  out.rhov = V.rho*V.v;
  out.rhoet = (V.p/(C.gamma-1.0)) + 0.5*V.rho*(V.u*V.u + V.v*V.v);
  return out;
}
//This function converts mach number into the primitive variables
primvar Mtoprim(double M, constants C)
{
  primvar answer;
  double psi = 1.0 + ((C.gamma - 1.0)/2.0) * M*M;//unitless
  double T = C.T0/psi;// K
  double R = 287.0;//J/(kg K)
  answer.p = ((C.p0*1000)/pow(psi,C.gamma/(C.gamma - 1.0)));//Pa
  answer.rho = answer.p/(R * T);
  answer.u = M * sqrt(C.gamma * R * T);
  return answer;
}
/*******************************************************************************/
/************************** END CONVERSIONS ************************************/
/*******************************************************************************/

/*******************************************************************************/
/************************** EXACT SOLUTION *************************************/
/*******************************************************************************/
prim_exact exactsol(double A_x, constants C)//constants is a datastruct defined in the header
{
  /********************************* Local Definitions *********************************/
  using namespace std;
  double A_t = C.A_t;//m^2
  double A_bar = A_x/A_t;//unitless
  double tol = C.tol;//unitless
  double phi;//unitless
  double F = 2.0;//unitless
  double dFdM;//unitless
  double Mnew;//unitless
  double dM=100000.0;//unitless
  //must compute the mach number and then with that one can compute the state variables
  double M;
  /**************************** Create the initial guess of Mach ***********************/
  if (C.cond == true) M=0.1;
  else M=5.0;
  /******************************* Newton Raphson Loop *********************************/
  while (abs(dM) > tol)
  {
    phi = (2.0/(C.gamma + 1.0)) * (1.0 + ((C.gamma - 1.0)/2.0) * M*M);//unitless
    F = pow(phi,(C.gamma + 1.0)/(C.gamma - 1.0)) - A_bar*A_bar * M*M;//unitless
    dFdM = 2.0 * M * (pow(phi,2.0/(C.gamma - 1.0)) - A_bar*A_bar);//unitless
    Mnew = M - (F/dFdM); //unitless 
    dM = (Mnew-M);//Used to check convergence
    M = Mnew;
  }
  //return the primary variables V once the Mach number is known.
  /***************************** Compute State Variables *******************************/

  prim_exact answer;
  primvar ans = Mtoprim(M, C);
  answer.rho = ans.rho;
  answer.u = ans.u;
  answer.p = ans.p;
  answer.M = M;
  return answer;//returns a prim_exact data struct defined in the header file.
}
/*******************************************************************************/
/************************ END EXACT SOLUTION ***********************************/
/*******************************************************************************/
