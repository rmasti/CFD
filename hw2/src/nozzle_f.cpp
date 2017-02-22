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

/****************************** RESIDUALS ***********************************/
void compute_residuals(vector<consvar> &Resarr, vector<double> &Res, vector<double> &Linfnorm, vector<double> &L1norm, vector<double> &L2norm, vector<consvar> const &Uold, vector<consvar> const &Unew)
{
  vector<double> sumresidual(3,0.0);
  vector<double> sumsqrresidual(3,0.0);

  for (int i=0; i<Uold.size()-1; i++)
  {
    sumresidual[0] += Resarr[i].rho;
    sumresidual[1] += Resarr[i].rhou;
    sumresidual[2] += Resarr[i].rhoet;

    sumsqrresidual[0] += Resarr[i].rho*Resarr[i].rho;
    sumsqrresidual[1] += Resarr[i].rhou*Resarr[i].rhou;
    sumsqrresidual[2] += Resarr[i].rhoet*Resarr[i].rhoet;

    if (Linfnorm[0] > Resarr[i].rho) Linfnorm[0] = Resarr[i].rho;
    if (Linfnorm[1] > Resarr[i].rhou) Linfnorm[1] = Resarr[i].rhou;
    if (Linfnorm[2] > Resarr[i].rhoet) Linfnorm[2] = Resarr[i].rhoet;

  }

  double invN = 1.0/Uold.size();
  L1norm[0] = sumresidual[0]*invN;
  L1norm[1] = sumresidual[1]*invN;
  L1norm[2] = sumresidual[2]*invN;

  L2norm[0] = sqrt(sumsqrresidual[0]*invN);
  L2norm[1] = sqrt(sumsqrresidual[1]*invN);
  L2norm[2] = sqrt(sumsqrresidual[2]*invN);

  Res[0] = L2norm[0];
  Res[1] = L2norm[1];
  Res[2] = L2norm[2];
  
}



/**************************** ITERATION STEP ********************************/

void iteration_step(vector<fluxes> &F, vector<consvar> &Uold, vector<consvar> &Unew, vector<primvar> &Vold, vector<primvar> &Vnew, vector<double> const &XCarr, vector<double> const &Xarr, vector<double> &Marr,vector<consvar> &Resarr, constants C) 
{

  vector<double> dtvec;
  double dtmin=1000.0; //localdt false then use dtmin global
  for(int i=0; i<Xarr.size()-1; i++) 
  {
    double dttemp = compute_timestep(Vold, i, C);
    dtvec.push_back(dttemp);
    if (dtmin > dttemp) dtmin = dttemp;
  }
  
  for(int i=1; i<XCarr.size(); i++)
  {
    double dt;
    // Compute time step for the sell i;
    if (localdt==true) dt = dtvec[i];
    else dt = dtmin;
    // Now that you have the time step create a function that computes the volume
    
    vector<double> ALR(2,0.0);
    double volume = compute_volume(Xarr,i, ALR);
    
    //cout << "i= " << i << " XCarr[i]= " <<  XCarr[i] <<" Uold[i].rho= " << Uold[i].rho << " F[i].rhou= " << F[i].rhou << endl;
    double S = Vold[i].p*dAdx(XCarr[i]);
    double inv_volume = 1.0/volume;
    // Take the step
   
    // continuity
    Unew[i].rho = Uold[i].rho - dt*inv_volume*(F[i].rhou*ALR[1] - F[i-1].rhou*ALR[0]);
    // x-mtm
    Unew[i].rhou = Uold[i].rhou + S*dx*dt*inv_volume - dt*inv_volume*(F[i].rhouu_and_p*ALR[1] - F[i-1].rhouu_and_p*ALR[0]);
    // y-mtm
    Unew[i].rhov = 0.0; //may return to this later
    // energy
    Unew[i].rhoet = Uold[i].rhoet - dt*inv_volume*(F[i].rhouht*ALR[1] - F[i-1].rhouht*ALR[0]);
    //cout << "flux cont = " << F[i].rhou << " " <<  F[i+1].rhou << endl;

    Vnew[i] = constoprim(Unew[i], C);
    Marr[i] = primtoM(Vnew[i], C);
/*
    Resarr[i].rho = dt*inv_volume*abs(Unew[i].rho - Uold[i].rho);
    Resarr[i].rhou = dt*inv_volume*abs(Unew[i].rhou - Uold[i].rhou);
    Resarr[i].rhoet = dt*inv_volume*abs(Unew[i].rhoet - Uold[i].rhoet);
*/
    Resarr[i].rho = dt*inv_volume*abs(Unew[i].rho - Uold[i].rho);
    //cout << " RHO: " << Resarr[i].rho << endl;
    Resarr[i].rhou = dt*inv_volume*abs(Unew[i].rhou - Uold[i].rhou);
    Resarr[i].rhoet = dt*inv_volume*abs(Unew[i].rhoet - Uold[i].rhoet);

  //  cout << "residual " << Unew[i].rho << endl;
  }
}

double compute_volume(vector<double> const &Xarr, int i, vector<double> &ALR)
{
  double xleft = Xarr[i];
  double xright = Xarr[i+1];

  double Aleft = A_x(xleft);
  double Aright = A_x(xright);
  double Aavg = 0.5*(Aleft+Aright);
  ALR[0] = Aleft;
  ALR[1] = Aright;

  return abs(xright-xleft)*Aavg;
}

double compute_timestep(vector<primvar> const &Vold, int i, constants C)
{
  double a = sqrt(Vold[i].p*C.gamma/Vold[i].rho);
  return cfl*dx/(abs(Vold[i].u) + a);
}

/********************************* FLUXES ***********************************/



void compute_fluxes(vector<fluxes> &F, vector<consvar> const &U,  constants C) 
{
  vector<consvar> U_avg; //vector of cons var on all the faces
  reconstruct_U(U_avg, U);
  //Now that we have the face avgd values of U we can compute F and also include some artificial disipation terms to the flux

  //nu and P are used for evaluating at values near i (i-2->i+2)...
  vector<double> P = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double epsilon2;
  double epsilon4;

  //loop and fill the traditional flux from U
  for(int i = 0 ; i < U_avg.size(); i++) F.push_back(fluxcalc(U_avg[i], C));

  //get the aritificial flux
 
  vector<fluxes> dvec;
 
  artificial_viscosity(dvec, U, C);

   //cout << F.size() << " "<< dvec.size() << endl;
  for(int i = 0 ; i < F.size(); i++)
  {
    F[i].rhou += dvec[i].rhou;
    F[i].rhouu_and_p += dvec[i].rhouu_and_p;
    F[i].rhouht += dvec[i].rhouht;
    //cout << " Fluxes at i " << F[i].rhou << " " << F[i].rhouu_and_p << " " << F[i].rhouht << endl;
  }
}
//This takes a given consvar vector and returns a flux vector
fluxes fluxcalc(consvar const &U, constants C)
{
  fluxes F;
  F.rhou = (U.rhou);
  F.rhouu_and_p = (0.5*(3.0 - C.gamma)*(U.rhou*U.rhou/U.rho))+(C.gamma - 1.0)*U.rhoet;
  F.rhouht = U.rhoet*(U.rhou/U.rho) + (U.rhou/U.rho)*(C.gamma - 1.0)*(U.rhoet - (0.5*U.rhou*U.rhou/U.rho));
  return F;
}

//This finds the average U values at all the faces
void reconstruct_U(vector<consvar> &U_avg, vector<consvar> const &U)
{
  // This will start at i = 3 which is the first cell of the domain of interest
  for(int i = num_ghost_cells-1; i<U.size()-num_ghost_cells; i++)
  {
    consvar Utemp;
    Utemp.rho = 0.5*(U[i+1].rho+U[i].rho);
    Utemp.rhou = 0.5*(U[i+1].rhou+U[i].rhou);
    Utemp.rhov = 0.5*(U[i+1].rhov+U[i].rhov);
    Utemp.rhoet = 0.5*(U[i+1].rhoet+U[i].rhoet);
    U_avg.push_back(Utemp);
    //cout << " Size of interface values "<< U_avg.size() << " " << U_avg.back().rhou << " " << U_avg.front().rhou << endl;
  }
  //cout << " SIZE OF " << U_avg.size() << endl;
}

/****************************************************************************/

/*************************** ARTIFICIAL VISCOSITY ***************************/
//return a dvec of all the fluxes on all interfaces
void artificial_viscosity(vector<fluxes> &dvec, vector<consvar> const &U, constants C)
{
  vector <double> P(6,10.0);
  vector <double> nu(4, 0.0);
 // need to loop over the interior points just like Uavg and Flux 
  for(int i = num_ghost_cells; i<=U.size()-num_ghost_cells; i++)
  {
    primvar V_2 = constoprim(U[i-2], C);
    primvar V_1 = constoprim(U[i-1], C);
    primvar V = constoprim(U[i], C);
    primvar V1 = constoprim(U[i+1], C);
    primvar V2 = constoprim(U[i+2], C);
    primvar V3 = constoprim(U[i+3], C);

    nu[0] = abs((V.p - 2.0*V_1.p - V_2.p)/(V.p + 2.0*V_1.p + V_2.p));
    nu[1] = abs((V1.p - 2.0*V.p - V_1.p)/(V1.p + 2.0*V.p + V_1.p));
    nu[2] = abs((V2.p - 2.0*V1.p - V.p)/(V2.p + 2.0*V1.p + V.p));
    nu[3] = abs((V3.p - 2.0*V2.p - V1.p)/(V3.p + 2.0*V2.p + V1.p));
    
    double numax = max(nu[0], nu[1]);
    numax = max(numax, nu[2]); numax = max(numax, nu[3]);
    double epsilon2 = kappa2*numax;
    double epsilon4 = kappa4*max(0.0, kappa4-epsilon2);


    double a_1 = sqrt(V_1.p*C.gamma/V_1.rho); 
    double a = sqrt(V.p*C.gamma/V.rho); 

    double lambda_1 = abs(V_1.u)+a_1;
    double lambda = abs(V.u)+a;
    double lambda_avg = 0.5*(lambda+lambda_1);

    //cout << i << endl;
    dvec.push_back(compute_dflux(epsilon2, epsilon4, lambda_avg, i, U)); 
  }
  //cout << " numax " << numax << endl;
}
//compute d flux
fluxes compute_dflux( double const &epsilon2, double const &epsilon4, double const &lambda, int const &i, vector<consvar> const &U)
{
  fluxes d;
  fluxes D1;
  D1.rhou = lambda*epsilon2*(U[i+1].rho - U[i].rho);
  D1.rhouu_and_p = lambda*epsilon2*(U[i+1].rhou - U[i].rhou);
  D1.rhouht = lambda*epsilon2*(U[i+1].rhoet - U[i].rhoet);

  fluxes D3;
  D3.rhou =lambda*epsilon4*(U[i+2].rho - 3.0*U[i+1].rho + 3.0*U[i].rho - U[i-1].rho);
  D3.rhouu_and_p =lambda*epsilon4*(U[i+2].rhou - 3.0*U[i+1].rhou + 3.0*U[i].rhou - U[i-1].rhou);
  D3.rhouht =lambda*epsilon4*(U[i+2].rhoet - 3.0*U[i+1].rhoet + 3.0*U[i].rhoet - U[i-1].rhoet);

  d.rhou = D1.rhou - D3.rhou;
  d.rhouu_and_p = D1.rhouu_and_p - D3.rhouu_and_p;
  d.rhouht = D1.rhouht - D3.rhouht;
  return d;
}

/****************************************************************************/

/*************************** BOUNDARY CONDITIONS ****************************/
//This function extrapolates interior to the ghost cells
void extrapolate_to_ghost( vector<consvar> &Uarr, constants C)
{
  int end = Uarr.size()-1;
  //if (end == Uarr.size()) exit(1);
  //cout << "Made it here" << endl;
  for (int i = 0; i < num_ghost_cells ; i++)
  {
    //extrapolate left values
    int ileft = num_ghost_cells-i;
    int iright = (end-num_ghost_cells)+i;

    //cout << "ileft = " << ileft-1 << " iright = " << iright+1 << endl;

    //left

    Uarr[ileft-1].rho = 2.0*Uarr[ileft].rho - Uarr[ileft+1].rho;
    Uarr[ileft-1].rhou = 2.0*Uarr[ileft].rhou - Uarr[ileft+1].rhou;
    Uarr[ileft-1].rhov = 2.0*Uarr[ileft].rhov - Uarr[ileft+1].rhov;
    Uarr[ileft-1].rhoet = 2.0*Uarr[ileft].rhoet - Uarr[ileft+1].rhoet;
    


    //right

    Uarr[iright+1].rho = 2.0*Uarr[iright].rho - Uarr[iright-1].rho;
    Uarr[iright+1].rhou = 2.0*Uarr[iright].rhou - Uarr[iright-1].rhou;
    Uarr[iright+1].rhov = 2.0*Uarr[iright].rhov - Uarr[iright-1].rhov;
    Uarr[iright+1].rhoet = 2.0*Uarr[iright].rhoet - Uarr[iright-1].rhoet;
  }
}

void set_boundary_cond( vector<consvar> &U, constants C) 
{
  // All three vectors have the right size for the ghost cells
  int end = (U.size())-1;
  double M_1 = primtoM(constoprim(U[1],C), C);
  double M_2 = primtoM(constoprim(U[2],C), C);
  double M_0;
  //inflow
  double in = 2.0*M_1 - M_2;

  if (in < 0.0) M_0 = 0.1;
  else M_0 = in;
 
  U[0] = primtocons(Mtoprim(M_0, C), C);

  //outflow
  if (C.outflow == true)
  {
    //Then the outflow is supersonic

    U[end].rho = 2.0*U[end-1].rho - U[end-2].rho;
    U[end].rhou = 2.0*U[end-1].rhou - U[end-2].rhou;
    U[end].rhov = 2.0*U[end-1].rhov - U[end-2].rhov;
    U[end].rhoet = 2.0*U[end-1].rhoet - U[end-2].rhoet;

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
}

/*****************************************************************************/



/********************************* WRITE OUT ********************************/
void write_out(FILE* & file, vector<double> const &Aarr, vector<double> const &XCarr, vector<consvar> const &U, constants C)
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
/*****************************************************************************/

/********************************* INITIALIZE ********************************/
//CHECKED
void initialize( vector<double> &M, vector<consvar> &U, constants C)
{
  //ADD LAYERS OF GHOST CELLS TO ALL THREE VECTORS
 // M.push_back(M_x(xmax+dx)); //add one cell to the end
  //M.insert(M.begin(), M_x(xmin+dx)); //add one to the beginning

  for(int i = 0; i < M.size(); i++) 
  { 
    primvar Vtemp = Mtoprim(M[i], C);
    U.push_back(primtocons(Vtemp, C));
  }
  //cout << " size " << V.size() << " " << U.size() << " " << M.size() << endl;
}
/*******************************************************************************/

/********************************* SET GEOMETRY ********************************/
//This function sets the x interface coordinates, and also the x center coordinates
void set_geometry(vector<double> &Aarr, vector<double> & xinterface , vector<double> &xcenter, vector<double> &Marr)
{

  double x_int_left = xmin-num_ghost_cells*dx;
  double x_int_right = xmax;
  //cout << dx << endl;
  for (int i=0; i < (N-1) + 2*num_ghost_cells ; i++)
  {
    double x_c = x_int_left + dx*(i+0.5);
    double x_i = x_int_left + dx*i;
    xcenter.push_back(x_c);
    xinterface.push_back(x_i);
    Aarr.push_back(A_x(x_c));
    Marr.push_back(M_x(x_c));
    //cout << xcenter.back() << " "<< xinterface.back() <<  endl;
  }

  xinterface.push_back(x_int_right+num_ghost_cells*dx); //add the last value
  //cout << " X_center " << xcenter.size() << " " << xinterface.size() << endl;
  //cout << xcenter.back() << " "<< xinterface.back() <<  endl;
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
  double slope = (1.9 - 0.1)/(xmax_dom - xmin_dom);
  double b = 0.1-slope*xmin_dom;
  if (x < xmin_dom) return slope*xmin_dom+b;
  if (x > xmax_dom ) return slope*xmax_dom+b;
  else return slope*x+b;
}
/****************************** END SET GEOMETRY *******************************/


/********************************* CONVERSIONS *********************************/
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
  out.rhoet = (V.p/(C.gamma-1.0)) + 0.5*V.rho*(V.u*V.u );
  return out;
}
//This function converts mach number into the primitive variables
primvar Mtoprim(double M, constants C)
{
  primvar answer;
  double psi = 1.0 + 0.5*(C.gamma - 1.0) * M*M;//unitless
  double T = C.T0/psi;// K
  answer.p = ((C.p0*1000)*pow(psi,(C.gamma - 1.0)/C.gamma));//Pa
  answer.rho = answer.p/(R * T);
  answer.u = M * sqrt(C.gamma * answer.p / answer.rho);
  answer.v = 0.0;
  return answer;
}

double primtoM(primvar V, constants C)
{
  double a = sqrt(V.p*C.gamma/V.rho); 
  return V.u/a;
}
/************************** END CONVERSIONS ************************************/

/************************** EXACT SOLUTION *************************************/
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
/************************ END EXACT SOLUTION ***********************************/
