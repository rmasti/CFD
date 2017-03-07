/*
 * Computational Fluid Dynamics
 * Written By: Robert Masti
 * 1/27/2017
 * This file contains all function declarations that will be called in hw1_m.cpp which is the 
 * main cpp file. This is the function file, which has its function prototypes stored in the 
 * header file.
 */
#include <iostream>
//#include <cmath>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <vector>

#include "nozzle.hpp" //structure templates and func prototypes

using namespace std;

/****************************** RESIDUALS ***********************************/

// This function computes the norms of the residuals for all three eqns
void compute_norms(
    vector<double> &Linfnorm,      //Output L-infinity norm vector with size==numeqns
    vector<double> &L1norm,        //Output L-1 norm  ||
    vector<double> &L2norm,        //Output L-2 norm  ||
    vector<consvar> const &Resarr  //Input Residual array with 3 elements 
    )
{
  //initialize vectors that will be used in summing the residuals
  vector<double> sumresidual(3,0.0);
  vector<double> sumsqrresidual(3,0.0);

  //loop over all interior cells 
  for (int i=num_ghost_cells; i<(Resarr.size()-num_ghost_cells); i++)
  {
    //sum residual for the L-1 norm   
    sumresidual[0] += Resarr[i].rho;
    sumresidual[1] += Resarr[i].rhou;
    sumresidual[2] += Resarr[i].rhoet;

    //sum the square residual for the L2 norm
    sumsqrresidual[0] += Resarr[i].rho*Resarr[i].rho;
    sumsqrresidual[1] += Resarr[i].rhou*Resarr[i].rhou;
    sumsqrresidual[2] += Resarr[i].rhoet*Resarr[i].rhoet;

    //Determine the Linf norm by comparing residual size 
    //NOTE: abs not needed here becase Resarr[:] >= 0.0
    if (Linfnorm[0] > Resarr[i].rho) Linfnorm[0] = Resarr[i].rho;
    if (Linfnorm[1] > Resarr[i].rhou) Linfnorm[1] = Resarr[i].rhou;
    if (Linfnorm[2] > Resarr[i].rhoet) Linfnorm[2] = Resarr[i].rhoet;
  }
  //finish computing the L1 & L2 norm respectively
  double invN = 1.0/Resarr.size();
  L1norm[0] = sumresidual[0]*invN;
  L1norm[1] = sumresidual[1]*invN;
  L1norm[2] = sumresidual[2]*invN;
  L2norm[0] = sqrt(sumsqrresidual[0]*invN);
  L2norm[1] = sqrt(sumsqrresidual[1]*invN);
  L2norm[2] = sqrt(sumsqrresidual[2]*invN);
}
/****************************************************************************/
/**************************** ITERATION STEP ********************************/
/****************************************************************************/
void iteration_step(
    vector<consvar> &Unew,         //Output (fill) Conserved Variable vector
    vector<consvar> &Resarr,       //Output (fill) residual array
    double &dt,
    vector<double> &Marr,          //Output (fill) Mach array
    vector<fluxes> const &F,       //Input fluxes (with artificial diss)
    vector<consvar> const  &Uold,  //Input the old consvar vector
    vector<double> const &XCarr,   //Input the center cell coords
    vector<double> const &Xarr,    //Input the interface coords
    constants const C              //Input constants for conversions, etc.
    )
{
  ////////////////////// DETERMINE TIME STEP /////////////////////
 
  //keep local dt determined from stability criteria
  vector<double> dtvec;
  double dtmin=1000.0; 
  // loop over the interior points of interest
  for(int i=num_ghost_cells; i < (Uold.size()-num_ghost_cells); i++) 
  {
    //get dt from compute_timestep funtion at cell i
    double dttemp = compute_timestep(Uold, i, C);
    dtvec.push_back(dttemp);
    //set dtmin for global time stepping
    if (dtmin > dttemp) dtmin = dttemp;
  }
 
  ////////////////////// MAIN SPATIAL LOOP ////////////////////////
  for(int i=num_ghost_cells; i < Uold.size() - num_ghost_cells ; i++) 
  {
    // this iterator will be set initially to zero for the fluxes that
    // are defined left of the center of the cell
    int i_interior = i-num_ghost_cells;
    // The dt used for time stepping
    // Compute time step for the sell i;
    if (localdt==true) dt = dtvec[i_interior];
    else dt = dtmin;
    /////// Compute the volume //////
    // ALR is the area on the left and right of the cell
    vector<double> ALR(2,0.0); 
    // return the volume and updated ALR
    double volume = compute_volume(ALR, Xarr,i);
    //cout << setprecision(14) << Xarr[i] << endl;
    /////// Compute Source term //////
    // NOTE there is two ways to determine S
    // Converting a const to prim allows for easy access to Pressure
    primvar Vold = constoprim(Uold[i], C);
    //cout << setprecision(14) << Vold.rho << endl;
    double S = Vold.p*dAdx(XCarr[i]); 
    //cout << setprecision(14) << S << endl;
    //double S = Vold.p*((ALR[1]-ALR[0])/dx);
    // Compute inverse volume as it is divided throughout
    double inv_volume = 1.0/volume;
    ///////// TAKE TIME STEP ////////
    // continuity
    //cout << setprecision(14) << dt << endl;
    //cout << setprecision(14) << volume << endl;
    Unew[i].rho = Uold[i].rho - dt*inv_volume*(F[i_interior+1].rhou*ALR[1] - F[i_interior].rhou*ALR[0]);
    // x-mtm
    Unew[i].rhou = Uold[i].rhou + S*dx*dt*inv_volume - dt*inv_volume*(F[i_interior+1].rhouu_and_p*ALR[1] - F[i_interior].rhouu_and_p*ALR[0]);
    // energy
    Unew[i].rhoet = Uold[i].rhoet - dt*inv_volume*(F[i_interior+1].rhouht*ALR[1] - F[i_interior].rhouht*ALR[0]);

    cout << setprecision(14) << S*dx + (F[i_interior+1].rhouu_and_p*ALR[1] - F[i_interior].rhouu_and_p*ALR[0]) << endl;

    // Fill the residual array
    // continuity
    Resarr[i].rho = dt*inv_volume*fabs(Unew[i].rho - Uold[i].rho);
    // x - momentum 
    Resarr[i].rhou = dt*inv_volume*fabs(Unew[i].rhou - Uold[i].rhou);
    // energy
    Resarr[i].rhoet = dt*inv_volume*fabs(Unew[i].rhoet - Uold[i].rhoet);
  }
}

double compute_timestep(vector<consvar> const &Uold, int i, constants const C)
{
  primvar Vold = constoprim(Uold[i], C);
  double a = sqrt(Vold.p*C.gamma/Vold.rho);
  double val = C.cfl*dx/(fabs(Vold.u) + a);
  return  val;
}

double compute_volume(vector<double> &ALR, vector<double> const &Xarr, int i)
{
  double xleft = Xarr[i];
  double xright = Xarr[i+1];
  double Aleft = A_x(xleft);
  double Aright = A_x(xright);
  double Aavg = 0.5*(Aleft+Aright);
  ALR[0] = Aleft;
  ALR[1] = Aright;
  //return fabs(xright-xleft)*Aavg;
  return dx*A_x(xleft+dx/2);
}

/****************************************************************************/
/*************************** FLUX IT HARDCORE  ******************************/
/****************************************************************************/
void compute_fluxes(
    vector<fluxes> &F,             //Output fluxes for all IOI (Interf of Inter)
    vector<consvar> const &U,      //Input consvar at all center cell values 
    constants const C              //Constants used for conversions
    ) 
{
  vector<consvar> U_avg; //vector of cons var on all the faces IOI

  // Reconstruct U values on to the faces
  reconstruct_U(U_avg, U);

  // loop and fill the traditional flux from U simple find F from U_avg
  for(int i = 0 ; i < U_avg.size(); i++) F.push_back(fluxcalc(U_avg[i], C));

  // get the aritificial flux
  vector<fluxes> dvec; //d-flux vector at all the same interface values
  artificial_viscosity(dvec, U, C); 

  ///// Combine d fluxes and F fluxes /////
  for(int i = 0 ; i < F.size(); i++)
  {

    //cout << setprecision(14) <<  F[i].rhou << endl;
    F[i].rhou = F[i].rhou + dvec[i].rhou;
    F[i].rhouu_and_p = F[i].rhouu_and_p + dvec[i].rhouu_and_p;
    F[i].rhouht = F[i].rhouht + dvec[i].rhouht;

    //cout << setprecision(14) <<  dvec[i].rhouu_and_p << endl;
    //cout << "flux cont rhou = " << F[i].rhou << " rhouu = " <<  F[i].rhouu_and_p << " rhouht = "<< F[i].rhouht << endl;
    if (isnan(F[i].rhou) || isnan(F[i].rhouu_and_p) || isnan(F[i].rhouht))
    {
      cout << "FLUX IS NAN" << endl;
      exit(1); 
    }
  }
}

//This takes a given consvar vector and returns a flux vector
fluxes fluxcalc(
    consvar const &U,              //Input avg U on faces
    constants const C              //Input C for conversions
    )
{
  // Initialize flux var
  fluxes F;
  // Fill Flux var
  // continuity (flux = mtm)
  F.rhou = U.rhou;
  // x-mtm (flux = enrgy)
  F.rhouu_and_p = ((3.0 - C.gamma)/2.0)*U.rhou*U.rhou/U.rho 
    + (C.gamma - 1.0)*U.rhoet;
  // energy (flux = enthalpy)
  F.rhouht = (U.rhoet*U.rhou)/U.rho + 
    (U.rhou/U.rho)*((C.gamma-1.0)*U.rhoet -
        ((C.gamma-1.0)/2.0)*U.rhou*U.rhou/U.rho);
  return F;
}

//This finds the average U values at all the faces
void reconstruct_U(
    vector<consvar> &U_avg,        //Output Uavg on all IOI             
    vector<consvar> const &U       //Input U on all cells including ghosts
    )
{
  // This will start at i = 2 which is the behind the first cell of the interface of interest
  for(int i = num_ghost_cells-1; i < U.size()-num_ghost_cells; i++)
  {
    // Compute avgs
    consvar Utemp;
    Utemp.rho = 0.5*(U[i+1].rho + U[i].rho);
    Utemp.rhou = 0.5*(U[i+1].rhou + U[i].rhou);
    Utemp.rhov = 0.5*(U[i+1].rhov + U[i].rhov);
    Utemp.rhoet = 0.5*(U[i+1].rhoet + U[i].rhoet);
    U_avg.push_back(Utemp);
  }
}

/****************************************************************************/
/*************************** ARTIFICIAL VISCOSITY ***************************/
/****************************************************************************/
//return a dvec of all the fluxes on all interfaces
void artificial_viscosity(
    vector<fluxes> &dvec,          //Output dflux vector at all interfaces of interest         
    vector<consvar> const &U,      //Input the cell values of the conserved variables
    constants const C              //Input constants for conversions
    )
{
  // Pressure's and nu's  for determining epsilon 2
  vector <double> P(6,10.0);
  vector <double> nu(4, 0.0);

  //////// Loop over the cells to get the interface values at the i+1/2 interfaces
  for(int i = num_ghost_cells-1; i < (U.size())-num_ghost_cells; i++)
  {
    // Convert the U to prim variable for easy access to pressure
    primvar V_2 = constoprim(U[i-2], C);
    primvar V_1 = constoprim(U[i-1], C);
    primvar V = constoprim(U[i], C);
    primvar V1 = constoprim(U[i+1], C);
    primvar V2 = constoprim(U[i+2], C);
    primvar V3 = constoprim(U[i+3], C);

    // Compute nu for the 4 cells 
    nu[0] = fabs((V.p - 2.0*V_1.p + V_2.p)/(V.p + 2.0*V_1.p + V_2.p));
    nu[1] = fabs((V1.p - 2.0*V.p + V_1.p)/(V1.p + 2.0*V.p + V_1.p));
    nu[2] = fabs((V2.p - 2.0*V1.p + V.p)/(V2.p + 2.0*V1.p + V.p));
    nu[3] = fabs((V3.p - 2.0*V2.p + V1.p)/(V3.p + 2.0*V2.p + V1.p));
   
    
    double numax =mymax(mymax(nu[0], nu[1]),mymax(nu[2], nu[3]));
    
    //cout<< setprecision(14) << "nu  " << nu[0] << " "<<  nu[1] << " "<<  nu[2] << " "<<  nu[3] << " "<<  endl; 

    // Find the epsilon values
    double epsilon2 = kappa2*numax;
    double epsilon4 = mymax(0.0, kappa4-epsilon2);

    // determine the lambda at i and i+1 to get i+1/2
    double a = sqrt(V.p*C.gamma/V.rho);  
    double a1 = sqrt(V1.p*C.gamma/V1.rho); 
    double lambda = fabs(V.u)+a;
    double lambda1 = fabs(V1.u)+a1;
    double lambda_avg = 0.5*(lambda+lambda1);

    //cout << setprecision(14) << "epsilon 2 = "<< epsilon2 << " epsilon4 =  " << epsilon4 << " lambda_avg =  " << lambda_avg <<  endl;
    //cout << setprecision(14) << " lambda_avg = " << lambda_avg << endl;

    //add the dflux to the dflux vector
        
    if (isnan(lambda_avg)) 
    {
      cout << "LAMBDA IS NAN" << endl;
      exit(1);
    }

    dvec.push_back(compute_dflux(epsilon2, epsilon4, lambda_avg, i, U)); 

    //cout << setprecision(14) << "dvec values rhou = " << dvec.back().rhou << " rhouu = " << dvec.back().rhouu_and_p << " rhouht = " << dvec.back().rhouht << endl;
  }

}
//compute d flux
fluxes compute_dflux( 
    double const &epsilon2,        //Input epsilon2
    double const &epsilon4,        //Input epsilon4
    double const &lambda,          //Input lambda  
    int const &i,                  //Input i value
    vector<consvar> const &U       //Input the Uconsvar vector
    )
{
  // Initialize the three fluxes
  fluxes d; //returned flux vector
  fluxes D1; 
  fluxes D3;

  // Compute D1
  D1.rhou = lambda*epsilon2*(U[i+1].rho - U[i].rho);
  D1.rhouu_and_p = lambda*epsilon2*(U[i+1].rhou - U[i].rhou);
  D1.rhouht = lambda*epsilon2*(U[i+1].rhoet - U[i].rhoet);

  // Compute D3
  D3.rhou =lambda*epsilon4*(U[i+2].rho - 3.0*U[i+1].rho + 3.0*U[i].rho - U[i-1].rho);
  D3.rhouu_and_p =lambda*epsilon4*(U[i+2].rhou - 3.0*U[i+1].rhou + 3.0*U[i].rhou - U[i-1].rhou);
  D3.rhouht =lambda*epsilon4*(U[i+2].rhoet - 3.0*U[i+1].rhoet + 3.0*U[i].rhoet - U[i-1].rhoet);

  // Combine the two
  d.rhou =  D3.rhou - D1.rhou;
  d.rhouu_and_p = D3.rhouu_and_p - D1.rhouu_and_p;
  d.rhouht = D3.rhouht - D1.rhouht;
  return d;
}

/***************************************************************************/
/*************************** BOUNDARY CONDITIONS ***************************/
/***************************************************************************/
//This function extrapolates interior to the ghost cells
void extrapolate_to_ghost( 
    vector<consvar> &Uarr,         //Output the ghost cell interped vals          
    constants const C              //Input Constants for conversions
    )
{
  int end = Uarr.size()-1;
  for (int i = 0; i < num_ghost_cells ; i++)
  {
    //extrapolate left values
    int ileft = num_ghost_cells-i;
    int iright = (end-num_ghost_cells)+i;
    //left
    Uarr[ileft-1].rho = (2.0*Uarr[ileft].rho - Uarr[ileft+1].rho);
    Uarr[ileft-1].rhou = (2.0*Uarr[ileft].rhou - Uarr[ileft+1].rhou);
    Uarr[ileft-1].rhov = (2.0*Uarr[ileft].rhov - Uarr[ileft+1].rhov);
    Uarr[ileft-1].rhoet =(2.0*Uarr[ileft].rhoet - Uarr[ileft+1].rhoet);
    //right
    Uarr[iright+1].rho = (2.0*Uarr[iright].rho - Uarr[iright-1].rho);
    Uarr[iright+1].rhou = (2.0*Uarr[iright].rhou - Uarr[iright-1].rhou);
    Uarr[iright+1].rhov = (2.0*Uarr[iright].rhov - Uarr[iright-1].rhov);
    Uarr[iright+1].rhoet = (2.0*Uarr[iright].rhoet - Uarr[iright-1].rhoet);
  }
}

void set_boundary_cond( 
    vector<consvar> &U,            //Output very end cell values      
    constants const C              //Input Constants for conversions
    )
{
  int end = (U.size())-1;

  // Extrapolate the mach number

  double M_1 = primtoM(constoprim(U[num_ghost_cells],C), C);
  double M_2 = primtoM(constoprim(U[num_ghost_cells+1],C), C);
  double M_0 = 0.0;
  for(int i=0; i < num_ghost_cells; i++)
  {
    int ib = (num_ghost_cells - 1) - i;
    M_0 = 2.0*M_1 - M_2;
    U[ib] = primtocons(Mtoprim(M_0, C), C);
    M_2 = M_1;
    M_1 = M_0;
  }

  //////// inflow ////////

  
  //////// outflow ////////
  if (C.outflow == true)
  {
    //Then the outflow is supersonic
    U[end].rho = (2.0*U[end-1].rho - U[end-2].rho);
    U[end].rhou = (2.0*U[end-1].rhou - U[end-2].rhou);
    U[end].rhov = (2.0*U[end-1].rhov - U[end-2].rhov);
    U[end].rhoet = (2.0*U[end-1].rhoet - U[end-2].rhoet);
    //assume ideal 
  }
  else
  {
    //Then the outflow is subsonic
    U[end].rho = 2.0*U[end-1].rho - U[end-2].rho;
    U[end].rhou = 2.0*U[end-1].rhou - U[end-2].rhou;
    U[end].rhov = 2.0*U[end-1].rhov - U[end-2].rhov;
    U[end].rhoet = 2.0*U[end-1].rhoet - U[end-2].rhoet;
    primvar Vend_1 = constoprim(U[end-1],C);
    primvar Vend = constoprim(U[end],C);
    Vend.p = 2.0*C.pb - Vend_1.p;
    U[end] = primtocons(Vend, C);
  }

  //cout << "BC: begin: rho = " << U[0].rho << " rhou = "<< U[0].rhou << " rhoet = " << U[0].rhoet << endl;

  //cout << "BC: end: rho = " << U[end].rho << " rhou = "<< U[end].rhou << " rhoet = " << U[end].rhoet << endl;
}

/***************************************************************************/
/******************************* SET GEOMETRY ******************************/
/***************************************************************************/

//This function sets the x interface coordinates, and also the x center coordinates
void set_geometry(
    vector<double> &Aarr,          //Outputs Area at cell centers             
    vector<double> &xinterface,    //Output interface coords
    vector<double> &xcenter,       //Output xcenter coords
    vector<double> &Marr)          //Output Mach number at centers
{

  double x_int_left = xmin-num_ghost_cells*dx;// left most x coord
  double x_int_right = xmax; //right most x coord

  // Create geometry
  for (int i=0; i < (N) + 2*num_ghost_cells ; i++)
  {
    double x_i = x_int_left + dx*i;
    double x_c = x_i+dx/2;
    xcenter.push_back(x_c);
    xinterface.push_back(x_i);
    Aarr.push_back(A_x(x_c));
    Marr.push_back(M_x(x_c));
    //cout << "x_c = "<< x_c << " Mach Number = " << M_x(x_c) << endl;
  }
  // There is always one more interface coord than cell center coord
  xinterface.push_back(x_int_right+num_ghost_cells*dx); //add the last value
}

void initialize( 
    vector<consvar> &U,            //Output the initial conserved vars       
    vector<double> const &M,       //Input Mach numbers
    constants const C              //Input C for conversions
    )
{
  //loop over mach values defined at centers
  for(int i = 0; i < M.size(); i++) U.push_back(primtocons(Mtoprim(M[i], C), C));
}

double A_x(double x)
{
  if (x < xmin_dom || x > xmax_dom) return 1.0;
  else return 0.2 + 0.4*(1.0 + sin(M_PI*(x-0.5)));//function that governs Area for nozzle
}
double dAdx(double x)
{
  if (x < xmin_dom || x > xmax_dom) return 0.0;
  else return 0.4*M_PI*cos(M_PI*(x-0.5));//derivative of func above
}
double M_x(double x)
{
  double slope = (1.9 - 0.1)/(xmax_dom - xmin_dom);
  double b = 0.1-slope*xmin_dom;
  double out; 
  if (x < xmin_dom) out = slope*xmin_dom+b;
  else if (x > xmax_dom ) out = slope*xmax_dom+b;
  else out = slope*x+b;
  return out;
}

/***************************************************************************/
/******************************* CONVERSIONS  ******************************/
/***************************************************************************/


double primtoM(primvar V, constants const C)
{
  double a = sqrt(V.p*C.gamma/V.rho); 
  return V.u/a;
}
primvar constoprim(
    consvar U, 
    constants const C)
{
  primvar out;
  out.rho = U.rho;
  double invrho = 1.0/U.rho;
  out.u = U.rhou/U.rho;
  out.v = U.rhov*invrho;
  out.p = (C.gamma-1.0)*(U.rhoet - 0.5*U.rhou*U.rhou/U.rho);
  return out;
}

consvar primtocons(
    primvar V, 
    constants const C)
{
  consvar out;
  out.rho = V.rho;
  out.rhou = V.rho*V.u;
  out.rhov = V.rho*V.v;
  out.rhoet = (V.p/(C.gamma-1.0)) + 0.5*V.rho*(V.u*V.u );
  return out;
}

//This function converts mach number into the primitive variables
primvar Mtoprim(double M, constants const C)
{
  primvar answer;
  double psi = 1.0 + 0.5*(C.gamma - 1.0) * M*M;//unitless
  double T = C.T0/psi;// K
  answer.p = ((C.p0*1000)/pow(psi,C.gamma/(C.gamma - 1.0)));//Pa
  answer.rho = answer.p/(R * T);
  answer.u = M * sqrt(C.gamma * answer.p / answer.rho);
  answer.v = 0.0;
  //cout << " primvar p = "<< answer.p<<" rho = "<<answer.rho<<" u = "<< answer.u << " psi= " << psi <<endl; 
  return answer;
}


/***************************************************************************/
/******************************** WRITE OUT ********************************/
/***************************************************************************/
void write_res(
    FILE* &file, 
    vector<int> const &n, 
    vector<double> const &t, 
    vector<double> const &Res1, 
    vector<double> const &Res2, 
    vector<double> const &Res3
    )
{
  for(int i = 0; i < n.size(); i++) 
  {
    fprintf(file, "%i %e %e %e %e\n", n[i], t[i], Res1[i], Res2[i], Res3[i]);
  }
}

void write_out(
    FILE* & file, 
    vector<double> const &Aarr, 
    vector<double> const &XCarr, 
    vector<consvar> const &U, 
    constants const C
    )
{
  int gc = num_ghost_cells;
  for(int i = (num_ghost_cells); i < XCarr.size()-num_ghost_cells ; i++) 
  {
    primvar V = constoprim(U[i], C);
    double M = primtoM(V, C);
    fprintf(file, "%e %e %e %e %e %e %e %e %e\n",XCarr[i], Aarr[i], V.rho, 
        V.u, V.p/1000.0, M, U[i].rho, U[i].rhou, U[i].rhoet);
  } 
}

/***************************************************************************/
/******************************** EXACT SOL ********************************/
/***************************************************************************/


prim_exact exactsol(double A_x, constants const C)//constants is a datastruct defined in the header
{
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
  /************* Create the initial guess of Mach ***********************/
  if (C.cond == true) M=0.1;
  else M=5.0;
  /************* Newton Raphson Loop *********************************/
  while (fabs(dM) > tol)
  {
    phi = (2.0/(C.gamma + 1.0)) * (1.0 + ((C.gamma - 1.0)/2.0) * M*M);//unitless
    F = pow(phi,(C.gamma + 1.0)/(C.gamma - 1.0)) - A_bar*A_bar * M*M;//unitless
    dFdM = 2.0 * M * (pow(phi,2.0/(C.gamma - 1.0)) - A_bar*A_bar);//unitless
    Mnew = M - (F/dFdM); //unitless 
    dM = (Mnew-M);//Used to check convergence
    M = Mnew;
  }
  //return the primary variables V once the Mach number is known.
  /************** Compute State Variables *******************************/

  prim_exact answer;
  primvar ans = Mtoprim(M, C);
  answer.rho = ans.rho;
  answer.u = ans.u;
  answer.p = ans.p;
  answer.M = M;
  return answer;//returns a prim_exact data struct defined in the header file.
}

/***************************************************************************/
/******************************** COUT COLL ********************************/
/***************************************************************************/
    //cout << "dt * inv_volume " << dt*inv_volume << endl;
    //cout << " X " << Xarr[i+1] << " "<<  Xarr[i]<<" " << XCarr[i] << endl;
    //cout << "cont rho " << dt*inv_volume*(F[i_interior+1].rhou*ALR[1] - F[i_interior].rhou*ALR[0]) << " " << Uold[i].rho <<endl;

    //cout << "cont rhou " << dt*inv_volume*(F[i_interior+1].rhouu_and_p*ALR[1] - F[i_interior].rhouu_and_p*ALR[0]) << " " << S*dx*dt*inv_volume << " "<< Uold[i].rhou <<endl;
/*
    double a = fabs((dt*inv_volume*(F[i_interior+1].rhou*ALR[1] - F[i_interior].rhou*ALR[0]))) ;
        
    double b =fabs((S*dx*dt*inv_volume - dt*inv_volume*(F[i_interior+1].rhouu_and_p*ALR[1] - F[i_interior].rhouu_and_p*ALR[0]))); 

    double c = fabs((dt*inv_volume*(F[i_interior+1].rhouht*ALR[1] - F[i_interior].rhouht*ALR[0])));
    a = fabs(Uold[i].rho - a)/fabs(Uold[i].rho);

    b = fabs(Uold[i].rhou - b)/fabs(Uold[i].rhou);

    c = fabs(Uold[i].rhoet - c)/fabs(Uold[i].rhoet);
*/
    //cout << "flux to old dif i "<< i << " rho " << a << " rhou " << b << " rhoet " << c<<endl;
    //cout << "cont rhoetet " << dt*inv_volume*(F[i_interior+1].rhouht*ALR[1] - F[i_interior].rhouht*ALR[0]) << " " << Uold[i].rhoet <<endl;
    //cout << i << " " << Unew[i].rho << " " << Unew[i].rhou << " " << Unew[i]. rhov << " " << Unew[i].rhoet << endl;

    //cout << "flux cont = " << F[i].rhou << " " <<  F[i].rhouu_and_p << F[i].rhouht << endl;
    /*cout << " Residual : " << abs(Unew[i].rho - Uold[i].rho) << " "<< 
      abs(Unew[i].rhou - Uold[i].rhou) << " " << 
      Unew[i].rhoet - Uold[i].rhoet << endl;
      */

  //cout << "Unew size "<< Unew.size() << endl;

    //cout << "nu  " << nu[0] << " "<<  nu[1] << " "<<  nu[2] << " "<<  nu[3] << " "<<  endl; 

    //cout << epsilon2 << " " << epsilon4 << " " << lambda_avg <<  endl;

    //cout << "dvec values " << dvec.back().rhou << " " << dvec.back().rhouu_and_p << " " << dvec.back().rhouht << endl;
  //cout << " dvecsize " << dvec.size() << endl;
