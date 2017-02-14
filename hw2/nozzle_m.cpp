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
#include "nozzle.h"
using namespace std;

int main()
{

  //Choose which data you wish to generate:
  //isentropic exact soln
  //isentropic fvm
  //shock fvm
  isentropic();
  return 0;
}

void isentropic()
{
  ofstream dataout; //create an output stream note iomanip lib must be included to change prec
  dataout.open ("exactsol.txt");// the file to open is data.txt in the dir in which it is run
  /********************************************************************************/
  /************************* Set Constants and D structs **************************/
  /********************************************************************************/
  constants consts; //constants is a struct declared in the header file
  double x[N]; //array of x values from -1 to 1 as defined in the problem
  double Aarr[N]; //array of areas that will be passed to the rootfinder funct.
  primvar prvparr[N]; //Another data struct defined in header.
  consts.p0 = 300.0; //kPa
  consts.T0 = 600.0; //K
  consts.A_t = A_x(0.0);//should be m^2
  consts.cond = true;
  consts.p0 = 300.0; //kPa
  consts.tol = 1.0e-8;
  consts.gamma = 1.4;
  /********************************************************************************/
  /***************************** Loop Over N Segments *****************************/
  /********************************************************************************/
  x[0]=-1.0; Aarr[0]=A_x(x[0]); prvparr[0] = exactsol(Aarr[0],consts);
  dataout  << "x (m) "<< " rho (kg)" << " u (m/s)" << " p (kPa)" << " M " <<endl;
  dataout << std::fixed << std::setprecision(14)
    << x[0] << " "
    << prvparr[0].rho << " "
    << prvparr[0].u << " "
    << prvparr[0].p/1000 << " "
    << prvparr[0].M << endl;//note the std::fixed, and std::setprecision are needed to output all sigfigs
  for (int i = 0; i<N+2;i++)
  {
    if (i<N/2)
    {
      consts.cond = true; //use a subsonic initial M guess
    }
    else
    {
      consts.cond = false; //use a supersonic initial M guess
    }
    x[i+1] = x[i] + 2.0/(N-1.0); //Evenly spaced x values
    Aarr[i+1] = A_x(x[i+1]); //The area is defined in hw1_f.cpp
    prvparr[i+1] = exactsol(Aarr[i+1],consts); //pass the datastruct and Area val to return all primary variables.
    dataout << std::fixed << std::setprecision(14)
      << x[i+1] << " "
      << prvparr[i+1].rho << " "
      << prvparr[i+1].u << " "
      << prvparr[i+1].p/1000 << " "
      << prvparr[i+1].M <<endl; //write out data
  }
  cout << "Done see exactsol.txt" << endl;
  dataout.close();
}
