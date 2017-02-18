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

#include "nozzle.hpp" //structure templates and func prototypes

using namespace std;

/********************************* FLUXES ***********************************/
void compute_fluxes(vector<fluxes> &F, vector<consvar> const &U, vector<primvar> const &V, constants C)  
{
  vector<consvar> U_avg; //vector of cons var on all the faces
  reconstruct_U(U_avg, U);
  //Now that we have the face avgd values of U we can compute F and also include some artificial disipation terms to the flux

  //nu and P are used for evaluating at values near i (i-2->i+2)...
  vector<double> nu = {0.0, 0.0, 0.0, 0.0};
  vector<double> P = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  vector<fluxes> dvec;
  //Loop over the middle points minus the first and last two so this takes usto 297 because there are 301 faces 
  for(int i=2; i < N-1; i++)
  {
    P[0] = V[i-2].p; P[1] = V[i-1].p; P[2] = V[i].p; 
    P[3] = V[i+1].p; P[4] = V[i+2].p; P[5] = V[i+3].p;
    nu[0] = abs((P[2] - 2.0*P[1] - P[0])/(P[2] + 2.0*P[1] + P[0]));
    nu[1] = abs((P[3] - 2.0*P[2] - P[1])/(P[3] + 2.0*P[2] + P[1]));
    nu[2] = abs((P[4] - 2.0*P[3] - P[2])/(P[4] + 2.0*P[3] + P[2]));
    nu[3] = abs((P[5] - 2.0*P[4] - P[3])/(P[5] + 2.0*P[4] + P[3]));
    double numax = max(nu[0], nu[1]);
    numax = max(numax, nu[2]); numax = max(numax, nu[3]);
    double epsilon2 = kappa2*numax;
    double epsilon4 = max(0.0, (kappa4-epsilon2));
    double lambda1 = abs(V[i].u) + sqrt(C.gamma * V[i].p / V[i].rho);
    double lambda2 = abs(V[i+1].u) + sqrt(C.gamma * V[i+1].p / V[i+1].rho);
    double lambda = 0.5*(lambda1 + lambda2);
    fluxes d;
    d.rhou = lambda*(epsilon2*(U[i+1].rho - U[i].rho) - epsilon4*(U[i+2].rho - 3.0*U[i+1].rho + 3.0*U[i].rho - U[i-1].rho));
    d.rhouu_and_p = lambda*(epsilon2*(U[i+1].rhou - U[i].rhou) - epsilon4*(U[i+2].rhou - 3.0*U[i+1].rhou + 3.0*U[i].rhou - U[i-1].rhou));
    d.rhouht = lambda*(epsilon2*(U[i+1].rhoet - U[i].rhoet) - epsilon4*(U[i+2].rhoet - 3.0*U[i+1].rhoet + 3.0*U[i].rhoet - U[i-1].rhoet));
    fluxes Ftemp;
    fluxcalc(Ftemp, U_avg[i], C);
    Ftemp.rhou += d.rhou;
    Ftemp.rhouu_and_p += d.rhouu_and_p;
    Ftemp.rhouht += d.rhouht;
    F.push_back(Ftemp); 
    dvec.push_back(d);//Need to save d so that we can interpolate to ghost cells and add to F in the end
  }
  
  //Interpolate d for the first and last two faces
  dvec.insert(dvec.begin(), {2.0*dvec[0].rhou - dvec[1].rhou, 
      2.0*dvec[0].rhouu_and_p - dvec[1].rhouu_and_p, 
      2.0*dvec[0].rhouht - dvec[1].rhouht});
  dvec.insert(dvec.begin(), {2.0*dvec[0].rhou - dvec[1].rhou, 
      2.0*dvec[0].rhouu_and_p - dvec[1].rhouu_and_p, 
      2.0*dvec[0].rhouht - dvec[1].rhouht});
  dvec.push_back({2.0*dvec[dvec.size()-1].rhou - dvec[dvec.size()-2].rhou, 
      2.0*dvec[dvec.size()-1].rhouu_and_p - dvec[dvec.size()-2].rhouu_and_p, 
      2.0*dvec[dvec.size()-1].rhouht - dvec[dvec.size()-2].rhouht});
  dvec.push_back({2.0*dvec[dvec.size()-1].rhou - dvec[dvec.size()-2].rhou, 
      2.0*dvec[dvec.size()-1].rhouu_and_p - dvec[dvec.size()-2].rhouu_and_p, 
      2.0*dvec[dvec.size()-1].rhouht - dvec[dvec.size()-2].rhouht});
  
  // Now take interpolated d's and add them to the flux determined from the average conserved var at faces for both the beginning and the end of vector F
  fluxes Ftemp;
  fluxcalc(Ftemp, U_avg[1], C);
  F.insert(F.begin(), {
       Ftemp.rhou + dvec[1].rhou,
       Ftemp.rhouu_and_p + dvec[1].rhouu_and_p,
       Ftemp.rhouht + dvec[1].rhouht
    });
  fluxcalc(Ftemp, U_avg[0], C);
  F.insert(F.begin(), {
       Ftemp.rhou + dvec[0].rhou,
       Ftemp.rhouu_and_p + dvec[0].rhouu_and_p,
       Ftemp.rhouht + dvec[0].rhouht
    });
  fluxcalc(Ftemp, U_avg[U_avg.size()-2], C);
  F.push_back({
       Ftemp.rhou + dvec[dvec.size()-2].rhou,
       Ftemp.rhouu_and_p + dvec[dvec.size()-2].rhouu_and_p,
       Ftemp.rhouht + dvec[dvec.size()-2].rhouht
    });
  fluxcalc(Ftemp, U_avg[U_avg.size()-1], C);
  F.push_back({
       Ftemp.rhou + dvec[dvec.size()-1].rhou,
       Ftemp.rhouu_and_p + dvec[dvec.size()-1].rhouu_and_p,
       Ftemp.rhouht + dvec[dvec.size()-1].rhouht
    });
  //interpolate at the end for the value of d on the left face
}

//This takes a given consvar vector and returns a flux vector
void fluxcalc(fluxes &F, consvar const &U, constants C)
{
  F.rhou = (U.rhou);
  F.rhouu_and_p = (0.5*(3.0 - C.gamma)*(U.rhou*U.rhou/U.rho));
  F.rhouht = (U.rhoet*(U.rhou/U.rho) + 
      (U.rhou/U.rho)*(C.gamma - 1.0)*(U.rhoet - 
          (0.5*U.rhou*U.rhou/U.rho)));
}

//This finds the average U values at all the faces
void reconstruct_U(vector<consvar> &U_avg, vector<consvar> const &U)
{
  for(int i = 0; i< U.size()-1; i++)
  {
    consvar Utemp;
    Utemp.rho = 0.5*(U[i+1].rho+U[i].rho);
    Utemp.rhou = 0.5*(U[i+1].rhou+U[i].rhou);
    Utemp.rhov = 0.5*(U[i+1].rhov+U[i].rhov);
    Utemp.rhoet = 0.5*(U[i+1].rhoet+U[i].rhoet);
    U_avg.push_back(Utemp);
  }
}

/****************************************************************************/

/*************************** BOUNDARY CONDITIONS ****************************/
void set_boundary_cond(vector<double> &M, vector<primvar> &V, vector<consvar> &U, constants C) 
{
  // All three vectors have the right size for the ghost cells
  int end = (M.size())-1;
  //inflow
  double in = 2.0*M[1] - M[2];

  if (in < 0.0) M[0] = 0.1;
  else M[0] = in;
  
  V[0] = Mtoprim(M[0], C);
  U[0] = primtocons(V[0], C);

  //outflow
  if (C.outflow == true)
  {
    //Then the outflow is supersonic
    V[end].rho = 2.0*V[end-1].rho - V[end-2].rho;
    V[end].u = 2.0*V[end-1].u - V[end-2].u;
    V[end].v = 2.0*V[end-1].v - V[end-2].v;
    V[end].p = 2.0*V[end-1].p - V[end-2].p;
    U[end] = primtocons(V[end], C);
    //assume ideal 
    double c = sqrt(C.gamma*V[end].p/V[end].rho);
    M[end] = V[end].u/c;
  }
  else
  {
    //Then the outflow is subsonic
    V[end].rho = 2.0*V[end-1].rho - V[end-2].rho;
    V[end].u = 2.0*V[end-1].u - V[end-2].u;
    V[end].v = 2.0*V[end-1].v - V[end-2].v;
    V[end].p = 2.0*C.pb - V[end-1].p;
    U[end] = primtocons(V[end], C);
    double c = sqrt(C.gamma*V[end].p/V[end].rho);
    M[end] = V[end].u/c; 
  }
}

/*****************************************************************************/



/********************************* WRITE OUT ********************************/
void write_out(FILE* & file, vector<double> const &Aarr, vector<double> const &Xarr, 
    vector<primvar> const &V, vector<double> const &M, vector<consvar> const &U)
{
  int gc = num_ghost_cells;
  for(int i = 0; i < Xarr.size() ; i++) 
  {
    fprintf(file, "%e %e %e %e %e %e %e %e %e\n",Xarr[i], Aarr[i], V[i+gc].rho, 
        V[i+gc].u, V[i+gc].p/1000.0, M[i+gc], U[i+gc].rho, U[i+gc].rhou, U[i+gc].rhoet);
  } 
}
/*****************************************************************************/

/********************************* INITIALIZE ********************************/
void initialize(vector<primvar> &V, vector<double> &M, vector<consvar> &U, constants C)
{

  //ADD LAYERS OF GHOST CELLS TO ALL THREE VECTORS
  M.push_back(0.0);
  M.insert(M.begin(), 0.0);

  for(int i = 0; i < M.size(); i++) 
  { 
    V.push_back(Mtoprim(M[i], C));
    U.push_back(primtocons(V[i], C));
  }

}
/*******************************************************************************/

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
  answer.p = ((C.p0*1000)/pow(psi,C.gamma/(C.gamma - 1.0)));//Pa
  answer.rho = answer.p/(R * T);
  answer.u = M * sqrt(C.gamma * R * T);
  answer.v = 0.0;
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
