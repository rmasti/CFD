/* 
 * Computational Fluid Dynamics
 * This code will compute the exact isentropic solution to the 1-D axisymmetric converging-diverging supersonic nozzle.
 * Written By: Robert Masti
 * 01/27/2017
 * This is the main file that will that will use a header file containing to function prototypes, structure declarations, symbolic constants, class declaration, etc.. And will use a seperate cpp file for the function declarations. 
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

#include "nozzle.hpp"
FILE *fp1;
FILE *fp2;
FILE *fp3;

using namespace std;

int main()
{

  //Choose which data you wish to generate:
  //isentropic exact soln
  //isentropic fvm
  //shock fvm
 
  constants consts; //constants is a struct declared in the header file
  consts.p0 = 300.0; //kPa
  consts.T0 = 600.0; //K
  consts.A_t = A_x(0.0);//should be m^2
  consts.cond = true;
  consts.tol = 1.0e-8;
  consts.gamma = 1.4;
  double dx = (xmax-xmin)/N;
  consts.nmax = 50000;
  consts.outflow = true;

  output_file_headers();

  isentropicExact(consts);
  isentropic(consts);
  fclose(fp1);
  cout << "Done see data directory" << endl;
  return 0;
}

void isentropic(constants C)
{
  vector<primvar> Varr;
  vector<consvar> Uarr;
  vector<double> Aarr;
  vector<double> Xarr;
  vector<double> dAdxarr;
  vector<double> Marr_i;

  //These two functions fill in the vectors that were defined above for both double vectors and also primvar struct defined in the header file
  set_geometry(Aarr, Xarr, dAdxarr, Marr_i);
  initialize(Varr, Marr_i , Uarr, C);

  set_boundary_cond(Marr_i, Varr, Uarr, C);

  write_out(fp2, Aarr, Xarr, Varr, Marr_i, Uarr);

  vector<primvar> Vold=Varr;
  vector<consvar> Uold=Uarr;
  vector<double> Res;
  Res.push_back(10.0); //rho
  Res.push_back(10.0);//rhou
  Res.push_back(10.0);//rhoet
  
  vector<fluxes> Farr;//defined on the cell faces
  compute_fluxes(Farr, Uarr, Varr, C);
 /* for (int n = 0; n < C.nmax; n++)
  {
  
  }*/
}


void isentropicExact(constants C)
{
  /************************* Set Constants and D structs **************************/
  double x[N+1]; //array of x values 
  double Aarr[N+1]; //array of areas that will be passed to the rootfinder funct.
  prim_exact prvparr[N+1]; //Another data struct defined in header.

  /***************************** Loop Over N Segments *****************************/
  for (int i = 0; i<N+1;i++)
  {
    if (i<N/2) C.cond = true; //use a subsonic initial M guess
    else C.cond = false; //use a supersonic initial M guess
    x[i] = xmin + i*(xmax-xmin)/(N); //Evenly spaced x values
    Aarr[i] = A_x(x[i]); //The area is defined in hw1_f.cpp
    prvparr[i] = exactsol(Aarr[i],C); //pass the datastruct and Area val to return all primary variables
    fprintf(fp1, "%f %e %e %e %f\n", x[i], prvparr[i].rho, prvparr[i].u, prvparr[i].p/1000.0, prvparr[i].M);
  }
}

/****************************** FILE HANDLING **********************************/
void output_file_headers()
{
  fp1 = fopen("data/exactsol.dat", "w");
  fprintf(fp1, "TITLE = \"Isentropic Exact Solution\"\n");
  fprintf(fp1, "variables=\"x(m)\"\"rho(kg/m^3)\"\"u(m/s)\"\"p(kPa)\"\"M()\"\n");

  fp2 = fopen("data/q1Dnozzle.dat","w");
  fprintf(fp2, "TITLE = \"Quasi-1D Nozzle Solution\"\n");
  fprintf(fp2, "variables=\"x(m)\"\"Area(m^2)\"\"rho(kg/m^3)\"\"u(m/s)\"\"Press(kPa)\"\"Mach\"\"U1\"\"U2\"\"U3\"\n");

  fp3 = fopen("data/q1history.dat","w");
  fprintf(fp3, "TITLE = \"Quasi-1D Residual History\"\n");
  fprintf(fp1,"variables=\"Iteration\"\"Time(s)\"\"Res1\"\"Res2\"\"Res3\"\n");
}
/***************************** END FILE HANDLING *******************************/

