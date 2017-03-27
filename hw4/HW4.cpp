// AOE6145 CFD HW4
// Euler equation for a quasi-1D nozzle
// Upwind Schemes ( 1st&2nd order of accuracy, Van Leer & Roe, Flux(slope) limiters )

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision

using namespace std;

const int flag_unitest = 1;     //Global flag; set 1 to turn unit test on, set 0 to turn it off

//===============================================================================================================//
//                                            Declare Functions                                                  //
//===============================================================================================================//
void BulidMesh(double DomainLeft, double DomainRight, int MeshSize, double *ptr_CellCenter, double *ptr_CellEdge);
double Calc_Area(double NozzleSpan);
double Calc_dAdx(double NozzleSpan);
void Ini_Mach(double NozzleSpan[], int MeshSize, double *ptr_Mach);
double Calc_M2phi(double Mach, double gamma);
double Calc_u2phi(double u, double T0, double gamma, double R);
double Calc_phi2P(double phi, double P0, double gamma);
double Calc_phi2rho(double phi, double P0, double T0, double R, double gamma);
void Ini_P(int MeshSize, double P0, double phi[], double gamma, double *ptr_P);
void Ini_rho(int MeshSize, double P0, double T0, double phi[], double R, double gamma, double *ptr_rho);
void Ini_u(int MeshSize, double Mach[], double P[], double rho[], double gamma, double *ptr_u);
void ExactCalc(int n, double A[], double Ac[], double xc[], double gamma, double R, double P0, double T0, \
	       double *ptr_M_exact, double *ptr_rho_exact, double *ptr_u_exact, double *ptr_P_exact);
double ExactCalc_InitialGuess(double NozzleSpan);
double ExactCalc_NewtonIteration(double IniGuess, double par1, double par2, double par3, \
				 int MaxIter, double eps);
double Newton_Function(double Mach, double Area, double ThroatArea, double gamma);
double Newton_FunctionDer(double Mach, double Area, double ThroatArea, double gamma);
void Switch_Prm2Consrv(int MeshSize, double gamma, double Prm1[], double Prm2[], double Prm3[], \
		       double *ptr_Consrv1, double *ptr_Consrv2, double *ptr_Consrv3);
void Switch_Consrv2Prm(int MeshSize, double gamma, double Consrv1[], double Consrv2[], double Consrv3[], \
		       double *ptr_Prm1, double *ptr_Prm2, double *ptr_Prm3);
void Switch_Consrv2Flux(int MeshSize, int num_ghost, double gamma, double Consrv1_tomb[], double Consrv2_tomb[], \
			double Consrv3_tomb[], double *ptr_Flux1, double *ptr_Flux2, double *ptr_Flux3);
void GoTomb(int MeshSize, int num_ghost, double realVar[], double *ptr_tombVar);
void OutTomb(int MeshSize, int num_ghost, double tombVar[], double *ptr_realVar);
void Calc_lambda_max(int MeshSize, int num_ghost, double gamma,		\
		     double rho_tomb[], double u_tomb[], double P_tomb[], double *ptr_lambda_max);
void Add_ArtiDissipation(int MeshSize, int num_ghost, double lambda_max[], double P_tomb[], \
			 double Consrv1_tomb[], double Consrv2_tomb[], double Consrv3_tomb[], \
			 double *ptr_Flux1, double *ptr_Flux2, double *ptr_Flux3);
void Calc_S(int MeshSize, double P[], double dAcdx[], double *ptr_S1, double *ptr_S2, double *ptr_S3);
double Max(double a, double b);
double SIGN(double a, double b);
void InflowBC(int MeshSize, int num_ghost, double u_tomb[], double P0, double T0, double gamma, double R, \
	      double *ptr_rho_tomb, double *ptr_P_tomb);
void OutflowBC(int flag_isen, int MeshSize, int num_ghost, double Pback, double *ptr_P_tomb);
void Calc_Res(int MeshSize, double Flux[], double Area[], double Source[], double dx[], double *ptr_Res);
double Calc_L2Norm(int MeshSize, double X[]);
void Calc_upwind_PrimLR(int MeshSize, int num_ghost, double upwind_eps, double upwind_kappa,
			    double prm_tomb[], double *ptr_prm_L, double *ptr_prm_R);
void Calc_upwind_Flux_VL(int MeshSize, double gamma, double R, double rho_L[], double rho_R[], double u_L[], \
			 double u_R[], double P_L[], double P_R[], double *ptr_F1, double *ptr_F2, double *ptr_F3);
void Calc_upwind_Flux_Roe(int MeshSize, double gamma, double R, double rho_L[], double rho_R[], double u_L[], \
			  double u_R[], double P_L[], double P_R[], double *ptr_F1, double *ptr_F2, double *ptr_F3);
//------------------------------------------- * I/O Functions * ---------------------------------------------//
void OutputSingleArray(string FileName, int ArrSize, double x[], double Arr[], int iteration);
void OutputSingleVector(string FileName, int VecSize, vector<double> Vec);
void CoutTitle(string TitleName, string LineStyle, string LengthStyle);

//===============================================================================================================//
//                                              Main Functions                                                   //
//===============================================================================================================//
int main()
{
  //***********************************************************************************************************//
  //                                           Setting Variables                                               //
  //***********************************************************************************************************//
  int flag_localTest;                 //Local Flag; set 1 to turn unit test on, set 0 to turn it off
  
  const int flag_isen = 1;      // set 1 to solve Isentropic case ( supersonic outflow )
                                // set 0 to solve for Normal shock case ( subsonic outflow )
  const  int flag_flux = 1;     // set 0 for Jameson Damping, set 1 for Van Leer, set 2 for Roe
  const double P0 = 300.0*1000.0;        // stagnation pressure in Pa
  const double T0 = 600.0;               // stagnation temperature in K
  const double Pback = 120.0*1000.0;     // ambient pressure
  const double gamma = 1.4;     
  const double R = 287.058;              // specific gas constatin in J/kg*K

  const int n = 256;                     // Number of cells. Make sure n is dividable by 6, so that there is a...
                                         //...face at the throat and the end of the two straight sections
  const double xmin = -1.5;              // Domain left boundary
  const double xmax = 1.5;               // Domain right boundary

  const int num_ghost = 3;               // ghost cells used at each boundary.

  const int WriteInterval = 1e4;
  
  double upwind_eps = 1;              // set 0 for 1st order. set 1 for 2nd or higher order
  double upwind_kappa = -1;           // -1 for fully-upwinded scheme; 0 for upwind-biased scheme; 1/3 for
                                      // 3rd order upwind- biased scheme (on uniform meshes); 1/2 for
                                      // Leonad's QUICK scheme; 1 for central difference scheme
  
  double CFL = 0.1;                  // CFL number
  int max_iter = 1e6;                // max number of iteration
  double Final_res = 1e-12;          // final residual
  double iter = 0;                   // iteration counter
 
  double xc[n], *ptr_xc = xc;   // array for cell center locations
  double x[n+1], *ptr_x = x;    // array for cell edge locations
  double dx[n];                 // array for dx. mesh is uniform in this case, but it can be non-uniform.
  double Min_dx;                // smallest dx in the domain
  double dt;                    // time step
  double A[n+1];                // array for the Area size at cell edges
  double Ac[n];                 // array for the Area size at cell centers
  double dAcdx[n];              // array for the first derivative of A, respect to x
  double Volumn[n];             // array for volumn of each cell
  double M[n], *ptr_M = M;               // array for the numerical solution of Mach number at cell centers
  double rho[n], *ptr_rho = rho;         // array for the numerical solution of Density at cell centers
  double u[n], *ptr_u = u;               // array for the numerical solution of Velocity at cell centers
  double P[n], *ptr_P = P;               // array for the numerical solution of Pressure at cell centers
  double M_exact[n], *ptr_M_exact = M_exact;       // array for the numerical solution of Mach # at cell centers
  double rho_exact[n], *ptr_rho_exact = rho_exact; // array for the exact solution of Density at cell centers
  double u_exact[n], *ptr_u_exact = u_exact;       // array for the exact solution of Velocity at cell centers
  double P_exact[n], *ptr_P_exact = P_exact;       // array for the exact solution of Pressure at cell centers
  double U1_exact[n], *ptr_U1_exact = U1_exact;    // array for the exact solution of Conserved 1 at cell centers
  double U2_exact[n], *ptr_U2_exact = U2_exact;    // array for the exact solution of Conserved 2 at cell centers
  double U3_exact[n], *ptr_U3_exact = U3_exact;    // array for the exact solution of Conserved 3 at cell centers

  double phi[n], *ptr_phi = phi;    // array for the phi at cell centers. (for primitive variables calc)  
  double lambda_max[n+1], *ptr_lambda_max = lambda_max;       // largest eigenvalue at each cell edge. 
  double Max_lambda;                // largest lambda in the domain.
  
  double U1[n], *ptr_U1 = U1;       // First conserved variable
  double U2[n], *ptr_U2 = U2;       // Second conserved variable
  double U3[n], *ptr_U3 = U3;       // Third conserved variable
  
  double F1[n+1], *ptr_F1 = F1;     // First Flux
  double F2[n+1], *ptr_F2 = F2;     // Second Flux
  double F3[n+1], *ptr_F3 = F3;     // Third Flux

  double S1[n], *ptr_S1 = S1;       // The source term for the Mass equation (which is 0 here)
  double S2[n], *ptr_S2 = S2;       // The source term for the Momentum equation
  double S3[n], *ptr_S3 = S3;       // The source term for the Energy equation (which is 0 here)
  
  double Res1[n], *ptr_Res1 = Res1; // Residual for Continuity equation
  double Res2[n], *ptr_Res2 = Res2; // Residual for Momentum equation
  double Res3[n], *ptr_Res3 = Res3; // Residual for Energy equation
  double L2norm1_initial, L2norm2_initial, L2norm3_initial;     // initial L2 norms of Residuals
  double L2norm1, L2norm2, L2norm3;                             // L2 norms of Residuals
  vector<double> L2norm1_rel;       // Relative L2 norm for Mass equation
  vector<double> L2norm2_rel;       // Relative L2 norm for Momentum equation
  vector<double> L2norm3_rel;       // Relative L2 norm for Energy equation

  //***********************************************************************************************************//
  //                                      Defined Upwind Scheme Variables                                      //
  //***********************************************************************************************************//
  double rho_L[n+1], *ptr_rho_L = rho_L;
  double rho_R[n+1], *ptr_rho_R = rho_R;
  double u_L[n+1], *ptr_u_L = u_L;
  double u_R[n+1], *ptr_u_R = u_R;
  double P_L[n+1], *ptr_P_L = P_L;
  double P_R[n+1], *ptr_P_R = P_R;
  
  //***********************************************************************************************************//
  //                      Defined Tomb Variables (contain both the ACTUAL and GHOST cells)                     //
  //***********************************************************************************************************//
  double rho_tomb[n+num_ghost*2], *ptr_rho_tomb = rho_tomb;
  double u_tomb[n+num_ghost*2], *ptr_u_tomb = u_tomb;
  double P_tomb[n+num_ghost*2], *ptr_P_tomb = P_tomb;

  double U1_tomb[n+num_ghost*2], *ptr_U1_tomb = U1_tomb;
  double U2_tomb[n+num_ghost*2], *ptr_U2_tomb = U2_tomb;
  double U3_tomb[n+num_ghost*2], *ptr_U3_tomb = U3_tomb;

  double xc_tomb[n+num_ghost*2], *ptr_xc_tomb = xc_tomb;
  
  //delete all previous output files in the "OutputFiles" folder
  system("exec rm ./OutputFiles/*");  
  //***********************************************************************************************************//
  //                                                 Meshing                                                   //
  //***********************************************************************************************************//
  BulidMesh(xmin, xmax, n, ptr_xc, ptr_x);       // Creat the Mesh of the domain;
  for (int i=0; i<n; i++)
    {
      dx[i] = x[i+1]-x[i];
      A[i] = Calc_Area(x[i]);
      Ac[i] = Calc_Area(xc[i]);
      dAcdx[i] = Calc_dAdx(xc[i]);
      Volumn[i] = dx[i]*Ac[i];
    }
  A[n] = Calc_Area(x[n]);
  Min_dx = *min_element(begin(dx),end(dx));     // The min dx in the domain

  //***********************************************************************************************************//
  //                                              Initialization                                               //
  //***********************************************************************************************************//
  // initialize Mach number and then calculate rho, u, P using Mach relation
  Ini_Mach(xc, n, ptr_M);   
  for (int i=0; i<n; i++) {
    phi[i] = Calc_M2phi(M[i], gamma);
  }
  Ini_P(n, P0, phi, gamma, ptr_P);
  Ini_rho(n, P0, T0, phi, R, gamma, ptr_rho);
  Ini_u(n, M, P, rho, gamma, ptr_u);


  // Set Boundary Conditions
  //    Convert real variable to tomb variable.
  //    compute ghost cells using linear extrapolation
  GoTomb(n, num_ghost, xc, ptr_xc_tomb);
  GoTomb(n, num_ghost, u, ptr_u_tomb);
  GoTomb(n, num_ghost, rho, ptr_rho_tomb);
  GoTomb(n, num_ghost, P, ptr_P_tomb);
  
  InflowBC(n, num_ghost, u_tomb, P0, T0, gamma, R, ptr_rho_tomb, ptr_P_tomb);
  OutflowBC(flag_isen, n, num_ghost, Pback, ptr_P_tomb);

  //---Switch Primitive variables to Conservative variables------//
  Switch_Prm2Consrv(n+2*num_ghost, gamma, rho_tomb, u_tomb, P_tomb, ptr_U1_tomb, ptr_U2_tomb, ptr_U3_tomb);
  
  //---------------Compute Flux---------------//
  Calc_lambda_max(n, num_ghost, gamma, rho_tomb, u_tomb, P_tomb, ptr_lambda_max ); // calc lambda_max...
                                                                                     // ...at every cell edge 
  Max_lambda = *max_element(begin(lambda_max)+num_ghost,end(lambda_max)-num_ghost);// The max lambda in the domain
  
  if (flag_flux == 0) { // central flux + jampson damping
    Switch_Consrv2Flux(n, num_ghost, gamma, U1_tomb, U2_tomb, U3_tomb, ptr_F1, ptr_F2, ptr_F3);
    Add_ArtiDissipation(n, num_ghost, lambda_max, P_tomb, U1_tomb, U2_tomb, U3_tomb, ptr_F1, ptr_F2, ptr_F3);
  }
  
  else { // upwind fluxes

    Calc_upwind_PrimLR(n, num_ghost, upwind_eps, upwind_kappa, rho_tomb, ptr_rho_L, ptr_rho_R);
    Calc_upwind_PrimLR(n, num_ghost, upwind_eps, upwind_kappa, u_tomb, ptr_u_L, ptr_u_R);
    Calc_upwind_PrimLR(n, num_ghost, upwind_eps, upwind_kappa, P_tomb, ptr_P_L, ptr_P_R);

  }
  if (flag_flux == 1)  // upwind: Van Leer
    Calc_upwind_Flux_VL(n, gamma, R, rho_L, rho_R, u_L, u_R, P_L, P_R, ptr_F1, ptr_F2, ptr_F3);
  if (flag_flux == 2)  // upwind: Roe
    Calc_upwind_Flux_Roe(n, gamma, R, rho_L, rho_R, u_L, u_R, P_L, P_R, ptr_F1, ptr_F2, ptr_F3);
  
  //---------------Compute Source-------------//
  Calc_S(n, P, dAcdx, ptr_S1, ptr_S2, ptr_S3);

  //----Convert tomb variables to real variables------//
  OutTomb(n, num_ghost, U1_tomb, ptr_U1);
  OutTomb(n, num_ghost, U2_tomb, ptr_U2);
  OutTomb(n, num_ghost, U3_tomb, ptr_U3);
  
  //------------Calculate Exact solution-------------------//
  ExactCalc(n, A, Ac, xc, gamma, R, P0, T0, ptr_M_exact, ptr_rho_exact, ptr_u_exact, ptr_P_exact); 
  Switch_Prm2Consrv(n, gamma, rho_exact, u_exact, P_exact, ptr_U1_exact, ptr_U2_exact, ptr_U3_exact);
  
  //------------Calculate Initial residual and norm-----------------//
  Calc_Res(n, F1, A, S1, dx, ptr_Res1);
  Calc_Res(n, F2, A, S2, dx, ptr_Res2);
  Calc_Res(n, F3, A, S3, dx, ptr_Res3);

  L2norm1_initial = Calc_L2Norm(n, Res1);
  L2norm2_initial = Calc_L2Norm(n, Res2);
  L2norm3_initial = Calc_L2Norm(n, Res3);

  L2norm1_rel.push_back(L2norm1_initial/L2norm1_initial);
  L2norm2_rel.push_back(L2norm2_initial/L2norm2_initial);
  L2norm3_rel.push_back(L2norm3_initial/L2norm3_initial);
  
  //------------Print out initial L2 norm of residuals----------//
  
  cout<<"\n**************Iteration = "<<iter<<" **************"<<endl;
  CoutTitle("Residuals", "=", "short");
  cout<<" * dt:                         "<<dt<<endl;
  cout<<" * Mass equation residual:     "<<L2norm1<<endl;
  cout<<" * Momentum equation residual: "<<L2norm2<<endl;
  cout<<" * Energy equation residual:   "<<L2norm3<<endl;
  
  
  //-------------------------------------------- * Output Files * ---------------------------------------------//
  flag_localTest =1;
  if (flag_localTest ==1)
    {
      OutputSingleArray("Area", n+1, x, A, iter);
      OutputSingleArray("lambda_max_initial", n+1, x, lambda_max, iter);
      OutputSingleArray("dx", n, xc, dx, iter);
      OutputSingleArray("Volumn", n, xc, Volumn, iter);
      OutputSingleArray("rho_initial", n, xc, rho, iter);
      OutputSingleArray("u_initial", n, xc, u, iter);
      OutputSingleArray("P_initial", n, xc, P, iter);

      OutputSingleArray("rho_L_initial", n+1, x, rho_L, iter);
      OutputSingleArray("u_L_initial", n+1, x, u_L, iter);
      OutputSingleArray("P_L_initial", n+1, x, P_L, iter);
      OutputSingleArray("rho_R_initial", n+1, x, rho_R, iter);
      OutputSingleArray("u_R_initial", n+1, x, u_R, iter);
      OutputSingleArray("P_R_initial", n+1, x, P_R, iter);
      
      OutputSingleArray("rho_tomb_initial", n+2*num_ghost, xc_tomb, rho_tomb, iter);
      OutputSingleArray("u_tomb_initial", n+2*num_ghost, xc_tomb, u_tomb, iter);
      OutputSingleArray("P_tomb_initial", n+2*num_ghost, xc_tomb, P_tomb, iter);
      OutputSingleArray("U1_tomb_initial", n+2*num_ghost, xc_tomb, U1_tomb, iter);
      OutputSingleArray("U2_tomb_initial", n+2*num_ghost, xc_tomb, U2_tomb, iter);
      OutputSingleArray("U3_tomb_initial", n+2*num_ghost, xc_tomb, U3_tomb, iter);
      OutputSingleArray("U1_initial", n, xc, U1, iter);
      OutputSingleArray("U2_initial", n, xc, U2, iter);
      OutputSingleArray("U3_initial", n, xc, U3, iter);
      OutputSingleArray("F1_initial", n+1, x, F1, iter);
      OutputSingleArray("F2_initial", n+1, x, F2, iter);
      OutputSingleArray("F3_initial", n+1, x, F3, iter);
      OutputSingleArray("S2_initial", n, xc, S2, iter);
      OutputSingleArray("Res1_initial", n, xc, Res1, iter);
      OutputSingleArray("Res2_initial", n, xc, Res2, iter);
      OutputSingleArray("Res3_initial", n, xc, Res3, iter);
      OutputSingleArray("rho_exact", n, xc, rho_exact, iter);
      OutputSingleArray("u_exact", n, xc, u_exact, iter);
      OutputSingleArray("P_exact", n, xc, P_exact, iter);
      OutputSingleArray("U1_exact", n, xc, U1_exact, iter);
      OutputSingleArray("U2_exact", n, xc, U2_exact, iter);
      OutputSingleArray("U3_exact", n, xc, U3_exact, iter);
    }

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //>>>>> Debug test----initial ghost cell
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  flag_localTest =0;
  if (flag_unitest == 1 && flag_localTest ==1)
    {
      cout<<"\n***************Unit Test---Ghost Cell***************"<<endl;
      cout<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<endl;
      cout << " 6 Ghost Cells of u "<<endl;
      for (int i=0; i<num_ghost; i++) {
	cout <<u_tomb[i]<<"; ";
      }
      cout<<" next is "<<u[0]<<"~"<<u[1]<<"~"<<u[2]<<"-----";
      
      for (int i=n+num_ghost; i<n+2*num_ghost; i++) {
	cout <<u_tomb[i]<<"; ";
      }
      cout << "\n Compara Interior cells of real and tomb u at begin, middle and the end: "<<endl;
      cout <<u[0]<<"---"<<u_tomb[num_ghost]<<"   ";
      cout <<u[n/2]<<"---"<<u_tomb[num_ghost+n/2]<<"   ";
      cout <<u[n-1]<<"---"<<u_tomb[num_ghost+n-1]<<"   "<<endl;    
    }

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //>>>>>> Debug test-----Initial condition
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  flag_localTest =0;
  if (flag_unitest == 1 && flag_localTest ==1)
    {
      cout<<"\n***************Unit Test---Initial condition***************"<<endl;
      cout<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<endl;
      cout<<"Initial Mach numbers at the beginning, the throat(2) and the end of the nozzle: \n" \
	  <<setprecision(14)<<M[0]<<"; "<<M[n/2-1]<<"; "<<M[n/2]<<"; "<<M[n-1]<<" ||"<<endl;
      cout<<"Initial u at the beginning, the throat(2) and the end of the nozzle: \n" \
	  <<setprecision(14)<<u[0]<<"; "<<u[n/2-1]<<"; "<<u[n/2]<<"; "<<u[n-1]<<" ||"<<endl;
      cout<<"Initial P at the beginning, the throat(2) and the end of the nozzle: \n" \
	  <<setprecision(14)<<P[0]<<"; "<<P[n/2-1]<<"; "<<P[n/2]<<"; "<<P[n-1]<<" ||"<<endl;
      cout<<"Initial rho at the beginning, the throat(2) and the end of the nozzle: \n" \
	  <<setprecision(14)<<rho[0]<<"; "<<rho[n/2-1]<<"; "<<rho[n/2]<<"; "<<rho[n-1]<<" ||"<<endl;
      cout<<"Initial lambda_max at the beginning, the throat and the end of the nozzle: \n" \
	  <<setprecision(14)<<lambda_max[0]<<"; "<<lambda_max[n/2]<<"; "<<lambda_max[n]<<"; "<<endl;
    }      
  
  //***********************************************************************************************************//
  //                                              Time Marching                                                //
  //***********************************************************************************************************//
  
  //---------------initialize dt---------------//
  dt = CFL*Min_dx/Max_lambda;
  
  while (iter < max_iter)
    {
      iter++;
      for (int i=0; i<n; i++)
	{
	  U1[i] = U1[i] - (dt/Volumn[i])*Res1[i];
	  U2[i] = U2[i] - (dt/Volumn[i])*Res2[i];
	  U3[i] = U3[i] - (dt/Volumn[i])*Res3[i];
	}
      Switch_Consrv2Prm(n, gamma, U1, U2, U3, ptr_rho, ptr_u, ptr_P);
      
      //----------Set Boundary Conditions---------//
      //    Convert real variable to tomb variable.
      //    compute ghost cells using linear extrapolation
      GoTomb(n, num_ghost, u, ptr_u_tomb);
      GoTomb(n, num_ghost, rho, ptr_rho_tomb);
      GoTomb(n, num_ghost, P, ptr_P_tomb);
      
      InflowBC(n, num_ghost, u_tomb, P0, T0, gamma, R, ptr_rho_tomb, ptr_P_tomb);
      OutflowBC(flag_isen, n, num_ghost, Pback, ptr_P_tomb);
      
      //---Switch Primitive variables to Conservative variables------//
      // update boundary conditions for Conservative variables
      Switch_Prm2Consrv(n+2*num_ghost, gamma, rho_tomb, u_tomb, P_tomb, ptr_U1_tomb, ptr_U2_tomb, ptr_U3_tomb);
      
      //---------------Compute Flux---------------//
      Calc_lambda_max(n, num_ghost, gamma, rho_tomb, u_tomb, P_tomb, ptr_lambda_max);  // calc lambda_max...
                                                                                       // ...at every cell edge 
      Max_lambda = *max_element(begin(lambda_max)+num_ghost,end(lambda_max)-num_ghost);

      if (flag_flux == 0) { // central flux + jampson damping
	Switch_Consrv2Flux(n, num_ghost, gamma, U1_tomb, U2_tomb, U3_tomb, ptr_F1, ptr_F2, ptr_F3);
	Add_ArtiDissipation(n, num_ghost, lambda_max, P_tomb, U1_tomb, U2_tomb, U3_tomb, ptr_F1, ptr_F2, ptr_F3);
      }

      else { // upwind fluxes

	Calc_upwind_PrimLR(n, num_ghost, upwind_eps, upwind_kappa, rho_tomb, ptr_rho_L, ptr_rho_R);
	Calc_upwind_PrimLR(n, num_ghost, upwind_eps, upwind_kappa, u_tomb, ptr_u_L, ptr_u_R);
	Calc_upwind_PrimLR(n, num_ghost, upwind_eps, upwind_kappa, P_tomb, ptr_P_L, ptr_P_R);
	
      }
      if (flag_flux == 1)  // upwind: Van Leer
	Calc_upwind_Flux_VL(n, gamma, R, rho_L, rho_R, u_L, u_R, P_L, P_R, ptr_F1, ptr_F2, ptr_F3);
      if (flag_flux == 2)  // upwind: Roe
	Calc_upwind_Flux_Roe(n, gamma, R, rho_L, rho_R, u_L, u_R, P_L, P_R, ptr_F1, ptr_F2, ptr_F3);
      
      //---------------Compute Source-------------//
      Calc_S(n, P, dAcdx, ptr_S1, ptr_S2, ptr_S3);

      //----Convert tomb variables to real variables------//
      OutTomb(n, num_ghost, U1_tomb, ptr_U1);
      OutTomb(n, num_ghost, U2_tomb, ptr_U2);
      OutTomb(n, num_ghost, U3_tomb, ptr_U3);

      //------------Calculate residual and norm-----------------//
      Calc_Res(n, F1, A, S1, dx, ptr_Res1);
      Calc_Res(n, F2, A, S2, dx, ptr_Res2);
      Calc_Res(n, F3, A, S3, dx, ptr_Res3);
    
      L2norm1 = Calc_L2Norm(n, Res1);
      L2norm2 = Calc_L2Norm(n, Res2);
      L2norm3 = Calc_L2Norm(n, Res3);

      L2norm1_rel.push_back(L2norm1/L2norm1_initial);
      L2norm2_rel.push_back(L2norm2/L2norm2_initial);
      L2norm3_rel.push_back(L2norm3/L2norm3_initial);
      
      //>>>>> Output files evert WriteInterval steps (With DEBUG feature for iter =1 <<<<<<//
      //------------------------------------------ * Output Files * -----------------------------------------//
      flag_localTest =1;
      if ( ((1.0*iter/WriteInterval)-int(1.0*iter/WriteInterval) ==0) ||\
	   (flag_unitest==1 && flag_localTest==1 && iter==1) )
	{
	  OutputSingleArray("rho_iter", n, xc, rho, iter);
	  OutputSingleArray("u_iter", n, xc, u, iter);
	  OutputSingleArray("P_iter", n, xc, P, iter);
	  OutputSingleArray("rho_tomb_iter", n+2*num_ghost, xc_tomb, rho_tomb, iter);
	  OutputSingleArray("u_tomb_iter", n+2*num_ghost, xc_tomb, u_tomb, iter);
	  OutputSingleArray("P_tomb_iter", n+2*num_ghost, xc_tomb, P_tomb, iter);
	  OutputSingleArray("U1_tomb_iter", n+2*num_ghost, xc_tomb, U1_tomb, iter);
	  OutputSingleArray("U2_tomb_iter", n+2*num_ghost, xc_tomb, U2_tomb, iter);
	  OutputSingleArray("U3_tomb_iter", n+2*num_ghost, xc_tomb, U3_tomb, iter);
	  OutputSingleArray("U1_iter", n, xc, U1, iter);
	  OutputSingleArray("U2_iter", n, xc, U2, iter);
	  OutputSingleArray("U3_iter", n, xc, U3, iter);
	  OutputSingleArray("F1_iter", n+1, x, F1, iter);
	  OutputSingleArray("F2_iter", n+1, x, F2, iter);
	  OutputSingleArray("F3_iter", n+1, x, F3, iter);
	  OutputSingleArray("S2_iter", n, xc, S2, iter);
	  OutputSingleArray("Res1_iter", n, xc, Res1, iter);
	  OutputSingleArray("Res2_iter", n, xc, Res2, iter);
	  OutputSingleArray("Res3_iter", n, xc, Res3, iter);
	}

      //-------Break if reach the final residual----------//
      if ((L2norm1_rel.back() < Final_res) && (L2norm2_rel.back() < Final_res) && (L2norm3_rel.back() < Final_res))
	{
	  OutputSingleArray("rho_final", n, xc, rho, 0);
	  OutputSingleArray("u_final", n, xc, u, 0);
	  OutputSingleArray("P_final", n, xc, P, 0);
	  OutputSingleArray("rho_tomb_final", n+2*num_ghost, xc_tomb, rho_tomb, 0);
	  OutputSingleArray("u_tomb_final", n+2*num_ghost, xc_tomb, u_tomb, 0);
	  OutputSingleArray("P_tomb_final", n+2*num_ghost, xc_tomb, P_tomb, 0);
	  OutputSingleArray("U1_tomb_final", n+2*num_ghost, xc_tomb, U1_tomb, 0);
	  OutputSingleArray("U2_tomb_final", n+2*num_ghost, xc_tomb, U2_tomb, 0);
	  OutputSingleArray("U3_tomb_final", n+2*num_ghost, xc_tomb, U3_tomb, 0);
	  OutputSingleArray("U1_final", n, xc, U1, 0);
	  OutputSingleArray("U2_final", n, xc, U2, 0);
	  OutputSingleArray("U3_final", n, xc, U3, 0);
	  OutputSingleArray("F1_final", n+1, x, F1, 0);
	  OutputSingleArray("F2_final", n+1, x, F2, 0);
	  OutputSingleArray("F3_final", n+1, x, F3, 0);
	  OutputSingleArray("S2_final", n, xc, S2, 0);
	  OutputSingleArray("Res1_final", n, xc, Res1, 0);
	  OutputSingleArray("Res2_final", n, xc, Res2, 0);
	  OutputSingleArray("Res3_final", n, xc, Res3, 0);
	  break;
	}

      //---------------update dt---------------//
      dt = CFL*Min_dx/Max_lambda;
      
      // Print out L2 norm of residuals
      if ((1.0*iter/WriteInterval)-int(1.0*iter/WriteInterval) ==0) {
	cout<<"\n**************Iteration = "<<iter<<" **************"<<endl;
	CoutTitle("Residuals", "=", "short");
	cout<<" * dt:                         "<<dt<<endl;
	cout<<" * Mass equation residual:     "<<L2norm1<<endl;
	cout<<" * Momentum equation residual: "<<L2norm2<<endl;
	cout<<" * Energy equation residual:   "<<L2norm3<<endl;
      }
    }

  //>>>>>>>>>>>> Output Relative L2 norms <<<<<<<<<<<//
  OutputSingleVector("L2norm1_rel", iter, L2norm1_rel);
  OutputSingleVector("L2norm2_rel", iter, L2norm2_rel);
  OutputSingleVector("L2norm3_rel", iter, L2norm3_rel);

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //>>>>>> Debug test----Final
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  flag_localTest =0;
  if (flag_unitest == 1 && flag_localTest ==1)
    {
      cout<<"\n********************Unit Test*********************"<<endl;
      cout<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<endl;
      cout<<"Cell Center locations at the first, middle and the last cell: \n"\
	  <<setprecision(14)<<xc[0]<<"; "<<xc[n/2]<<"; "<<xc[n-1]<<" ||"<<endl;
      cout<<"Cell edge locations at the first, middle and the last cell: \n"\
	  <<setprecision(14)<<x[0]<<"; "<<x[n/2]<<"; "<<x[n]<<" ||"<<endl;
      cout<<"Nozzle area sizes at x = -1.5, -1, 0, 1, 1.5: \n"		\
	  <<setprecision(14)<<A[0]<<"; "<<A[n/6]<<"; "<<A[n/2]<<"; "<<A[5*n/6]<<"; "<<A[n]<<" ||"<<endl;
      cout<<"M_exact at the beginning, the throat(2) and the end of the nozzle: \n"\
	  <<setprecision(14)<<M_exact[0]<<"; "<<M_exact[n/2-1]<<"; "<<M_exact[n/2]<<"; "<<M_exact[n-1]<<" ||"<<endl;
    }

  //***********************************************************************************************************//
  //                                           Free Computer Memory                                            //
  //***********************************************************************************************************//
  ptr_xc = new double;               delete ptr_xc;               ptr_xc = NULL;
  ptr_x = new double;                delete ptr_x;                ptr_x = NULL;

  ptr_M = new double;                delete ptr_M;                ptr_M = NULL;
  ptr_rho = new double;              delete ptr_rho;              ptr_rho = NULL;
  ptr_u = new double;                delete ptr_u;                ptr_u = NULL;
  ptr_P = new double;                delete ptr_P;                ptr_P = NULL;

  ptr_M_exact = new double;          delete ptr_M_exact;          ptr_M_exact = NULL;
  ptr_rho_exact = new double;        delete ptr_rho_exact;        ptr_rho_exact = NULL;
  ptr_u_exact = new double;          delete ptr_u_exact;          ptr_u_exact = NULL;
  ptr_P_exact = new double;          delete ptr_P_exact;          ptr_P_exact = NULL;
  ptr_U1_exact = new double;         delete ptr_U1_exact;         ptr_U1_exact = NULL;
  ptr_U2_exact = new double;         delete ptr_U2_exact;         ptr_U2_exact = NULL;
  ptr_U3_exact = new double;         delete ptr_U3_exact;         ptr_U3_exact = NULL;
  
  ptr_phi = new double;              delete ptr_phi;              ptr_phi = NULL;
  
  ptr_U1 = new double;               delete ptr_U1;               ptr_U1 = NULL;
  ptr_U2 = new double;               delete ptr_U2;               ptr_U2 = NULL;
  ptr_U3 = new double;               delete ptr_U3;               ptr_U3 = NULL;

  ptr_F1 = new double;               delete ptr_F1;               ptr_F1 = NULL;
  ptr_F2 = new double;               delete ptr_F2;               ptr_F2 = NULL;
  ptr_F3 = new double;               delete ptr_F3;               ptr_F3 = NULL;

  ptr_S1 = new double;               delete ptr_S1;               ptr_S1 = NULL;
  ptr_S2 = new double;               delete ptr_S2;               ptr_S2 = NULL;
  ptr_S3 = new double;               delete ptr_S3;               ptr_S3 = NULL;
  
  ptr_Res1 = new double;             delete ptr_Res1;             ptr_Res1 = NULL;
  ptr_Res2 = new double;             delete ptr_Res2;             ptr_Res2 = NULL;
  ptr_Res3 = new double;             delete ptr_Res3;             ptr_Res3 = NULL;

  ptr_rho_L = new double;            delete ptr_rho_L;            ptr_rho_L = NULL;
  ptr_rho_R = new double;            delete ptr_rho_R;            ptr_rho_R = NULL;

  ptr_u_L = new double;              delete ptr_u_L;              ptr_u_L = NULL;
  ptr_u_R = new double;              delete ptr_u_R;              ptr_u_R = NULL;

  ptr_P_L = new double;              delete ptr_P_L;              ptr_P_L = NULL;
  ptr_P_R = new double;              delete ptr_P_R;              ptr_P_R = NULL;
  
  ptr_rho_tomb = new double;         delete ptr_rho_tomb;         ptr_rho_tomb = NULL;
  ptr_u_tomb = new double;           delete ptr_u_tomb;           ptr_u_tomb = NULL;
  ptr_P_tomb = new double;           delete ptr_P_tomb;           ptr_P_tomb = NULL;

  ptr_U1_tomb = new double;          delete ptr_U1_tomb;          ptr_U1_tomb = NULL;
  ptr_U2_tomb = new double;          delete ptr_U2_tomb;          ptr_U2_tomb = NULL;
  ptr_U3_tomb = new double;          delete ptr_U3_tomb;          ptr_U3_tomb = NULL;
  
  ptr_xc_tomb = new double;          delete ptr_xc_tomb;          ptr_xc_tomb = NULL;
}

//===============================================================================================================//
//                                          Other Numerical Functions                                            //
//===============================================================================================================//
void BulidMesh(double DomainLeft, double DomainRight, int MeshSize, double *ptr_CellCenter, double *ptr_CellEdge)
{
  double dx = (DomainRight - DomainLeft)/MeshSize; 
  for (int i = 0; i < MeshSize; i++)
    {
      *(ptr_CellCenter+i) = DomainLeft + (0.5+i)*dx;
      *(ptr_CellEdge+i) = DomainLeft + i*dx;
    }
  *(ptr_CellEdge+MeshSize) = DomainRight;  
}

double Calc_Area(double NozzleSpan)
{
  double pi = acos(-1);
  double Area;
  if (NozzleSpan < -1.0 || NozzleSpan > 1.0)
    Area = 1.0;
  if (NozzleSpan >= -1.0 && NozzleSpan <= 1.0) 
    Area = 0.2 + 0.4*(1.0+sin(pi*(NozzleSpan-0.5)));  
  return Area;
}

double Calc_dAdx(double NozzleSpan)
{
  double pi = acos(-1);
  double dAdx;
  if (NozzleSpan <= -1.0 || NozzleSpan >= 1.0)
    dAdx = 0.0;
  if (NozzleSpan > -1.0 && NozzleSpan < 1.0)
    dAdx = 0.4*pi*cos(pi*(NozzleSpan-0.5));
  return dAdx;
}

void Ini_Mach(double NozzleSpan[], int MeshSize, double *ptr_Mach)
{
  for (int i = 0; i < MeshSize; i++)
    {
      if (NozzleSpan[i] < -1.0)
	*(ptr_Mach+i) = 0.2;
      if (NozzleSpan[i] >= -1.0 && NozzleSpan[i] <= 1.0)
	*(ptr_Mach+i) = 0.9*NozzleSpan[i] + 1.0;
      if (NozzleSpan[i] > 1.0)
	*(ptr_Mach+i) = 1.8;
    }
}

double Calc_M2phi(double Mach, double gamma)
{
  return 1.0+(gamma-1.0)*Mach*Mach/2.0;
}

double Calc_u2phi(double u, double T0, double gamma, double R)
{
  return T0/( T0 - (gamma-1)*u*u/(2.0*gamma*R) );    //  phi = T0/T
}

double Calc_phi2P(double phi, double P0, double gamma)
{
  return P0/pow(phi, gamma/(gamma-1) );
}

double Calc_phi2rho(double phi, double P0, double T0, double R, double gamma)
{
  return P0/pow(phi, 1.0/(gamma-1) )/(R*T0);
}
  
void Ini_P(int MeshSize, double P0, double phi[], double gamma, double *ptr_P)
{
  for (int i = 0; i < MeshSize; i++) {
    *(ptr_P+i) = P0/pow(phi[i],gamma/(gamma-1));
  }
}

void Ini_rho(int MeshSize, double P0, double T0, double phi[], double R, double gamma, double *ptr_rho)
{
  double rho0 = P0/R/T0;
  for (int i = 0; i < MeshSize; i++) {
    *(ptr_rho+i) =  rho0/pow(phi[i],1.0/(gamma-1));
  }
}

void Ini_u(int MeshSize, double Mach[], double P[], double rho[], double gamma, double *ptr_u)
{
  for (int i &= 0; i < MeshSize; i++) {
    *(ptr_u+i) = Mach[i]*sqrt(gamma*P[i]/rho[i]);
  }
}
  
void ExactCalc(int n, double A[], double Ac[], double xc[], double gamma, double R, double P0, double T0, \
	       double *ptr_M_exact, double *ptr_rho_exact, double *ptr_u_exact, double *ptr_P_exact)
{
  double M0;             // Initial guess for Mach Number 
  double A0 = A[n/2];    // Area size for the throat section  
  double phi;            // Used to calculate the primitive variables  

  int MaxIter = 1000;      // Max number of iterations for Newton Iteration
  double eps = 1e-12;         // Convergence criteria for Newton Iteration
  for (int i = 0; i < n; i++)
    {
      M0 = ExactCalc_InitialGuess(xc[i]);   // Calculate the initial guess of Mach number at each cell center
      *(ptr_M_exact+i) = ExactCalc_NewtonIteration(M0,gamma,Ac[i],A0,MaxIter,eps);     
                // Newton Iteration to calculate Mach number for the given span.

      phi = 1+(gamma-1)**(ptr_M_exact+i)**(ptr_M_exact+i)/2;
      *(ptr_P_exact+i) = P0/pow(phi,gamma/(gamma-1));
      *(ptr_rho_exact+i) = *(ptr_P_exact+i)*phi/(R*T0);
      *(ptr_u_exact+i) = *(ptr_M_exact+i)*sqrt(gamma*R*T0/phi);
    }
}
  
double ExactCalc_InitialGuess(double NozzleSpan)
{
  double Mach;
  if (NozzleSpan <= 0.0)
    Mach = 0.5;
  if (NozzleSpan > 0.0) 
    Mach = 3;
  return Mach;
}

double ExactCalc_NewtonIteration(double IniGuess, double par1, double par2, double par3, int MaxIter, double eps)
// par1 = gamma; par2 = A; par3 = A0;
{
  double F;     // The function to be solved in the Newton Iteration
  double dFdM;  // The derivative of this function respect to Mach number
  double Mach;  // Mach number
  double dM;    // delta M
  
  F = Newton_Function(IniGuess, par2, par3, par1);
  dFdM = Newton_FunctionDer(IniGuess, par2, par3, par1);
  
  Mach = IniGuess;
  int iter = 0;        // Iteration counter
  while (abs(F) > eps)
    {
      dM = -F/dFdM;
      Mach = Mach + dM;
      F = Newton_Function(Mach, par2, par3, par1);
      dFdM = Newton_FunctionDer(Mach, par2, par3, par1);
      iter++;
      if (iter == MaxIter) 
	break;  
    }
  return Mach;
}

double Newton_Function(double Mach, double Area, double ThroatArea, double gamma)
{
  return pow((1.0+(gamma-1.0)*Mach*Mach/2.0)*2.0/(gamma+1.0), (gamma+1.0)/(gamma-1.0)) \
    - pow(Mach*Area/ThroatArea ,2.0);
}

double Newton_FunctionDer(double Mach, double Area, double ThroatArea, double gamma)
{
  return 2.0*Mach*(pow((1.0+(gamma-1.0)*Mach*Mach/2.0)*2.0/(gamma+1.0), 2.0/(gamma-1.0)) \
		   - pow(Area/ThroatArea ,2.0));
}

void Switch_Prm2Consrv(int MeshSize, double gamma, double Prm1[], double Prm2[], double Prm3[], \
		       double *ptr_Consrv1, double *ptr_Consrv2, double *ptr_Consrv3)
{ // Prm1 = rho, Prm2 = u, Prm3 = P;    Consrv1 = rho, Consrv2 = rho*u, Consrv3 = rho*et;
  for (int i = 0; i < MeshSize; i++)
    {
      *(ptr_Consrv1+i) = Prm1[i];
      *(ptr_Consrv2+i) = Prm1[i]*Prm2[i];
      *(ptr_Consrv3+i) = Prm3[i]/(gamma-1) + Prm1[i]*Prm2[i]*Prm2[i]/2;
    } 
}

void Switch_Consrv2Prm(int MeshSize, double gamma, double Consrv1[], double Consrv2[], double Consrv3[], \
		       double *ptr_Prm1, double *ptr_Prm2, double *ptr_Prm3)
{ // Prm1 = rho, Prm2 = u, Prm3 = P;    Consrv1 = rho, Consrv2 = rho*u, Consrv3 = rho*et;
  for (int i = 0; i < MeshSize; i++)
    {
      *(ptr_Prm1+i) = Consrv1[i];
      *(ptr_Prm2+i) = Consrv2[i]/Consrv1[i];
      *(ptr_Prm3+i) = (gamma-1)*(Consrv3[i] - Consrv2[i]*Consrv2[i]/(Consrv1[i]*2.0));
    }
}

void Switch_Consrv2Flux(int MeshSize, int num_ghost, double gamma, double Consrv1_tomb[], double Consrv2_tomb[], \
			double Consrv3_tomb[], double *ptr_Flux1, double *ptr_Flux2, double *ptr_Flux3)
{
  double Uc1[MeshSize+1], Uc2[MeshSize+1], Uc3[MeshSize+1]; 
  for (int i=0; i < MeshSize+1; i++)
    {
      Uc1[i] = 0.5*(Consrv1_tomb[num_ghost+i-1]+Consrv1_tomb[num_ghost+i]);
      Uc2[i] = 0.5*(Consrv2_tomb[num_ghost+i-1]+Consrv2_tomb[num_ghost+i]);
      Uc3[i] = 0.5*(Consrv3_tomb[num_ghost+i-1]+Consrv3_tomb[num_ghost+i]);
    }
  for (int i =0; i < MeshSize+1; i++)
    {
      *(ptr_Flux1+i) = Uc2[i];
      *(ptr_Flux2+i) = 0.5*(3-gamma)*Uc2[i]*Uc2[i]/Uc1[i] + (gamma-1)*Uc3[i];
      *(ptr_Flux3+i) = gamma*Uc2[i]*Uc3[i]/Uc1[i] - 0.5*(gamma-1)*Uc2[i]*Uc2[i]*Uc2[i]/(Uc1[i]*Uc1[i]);
      //*(ptr_Flux3+i) = Uc3[i]*Uc2[i]/Uc1[i]+(Uc2[i]/Uc1[i])*((gamma-1)*Uc3[i] - 0.5*(gamma-1)*Uc2[i]*Uc2[i]/Uc1[i]);
    }
}

void GoTomb(int MeshSize, int num_ghost, double realVar[], double *ptr_tombVar)
{
  double *ghost_Lcursor;  //point address for a left ghost cell to be manipulated
  double *ghost_Rcursor;  //point address for a right ghost cell to be manipulated

  copy(realVar, realVar+MeshSize, ptr_tombVar+num_ghost);

  for (int i = 0; i < num_ghost; i++)
    {
      ghost_Lcursor = ptr_tombVar + num_ghost - i - 1;       // Start from the left closest ghost cell
      ghost_Rcursor = ptr_tombVar + num_ghost + MeshSize + i;// Start from the right closest ghost cell
      *(ghost_Lcursor) = 2.0**(ghost_Lcursor+1) - *(ghost_Lcursor+2);  
      *(ghost_Rcursor) = 2.0**(ghost_Rcursor-1) - *(ghost_Rcursor-2);
    }
  ghost_Lcursor = new double;   delete ghost_Lcursor;   ghost_Lcursor = NULL;
  ghost_Rcursor = new double;   delete ghost_Rcursor;   ghost_Rcursor = NULL;
}

void OutTomb(int MeshSize, int num_ghost, double tombVar[], double *ptr_realVar)
{
  copy(tombVar+num_ghost, tombVar+num_ghost+MeshSize, ptr_realVar);
}
  
void Calc_lambda_max(int MeshSize, int num_ghost, double gamma, \
		     double rho_tomb[], double u_tomb[], double P_tomb[], double *ptr_lambda_max)
{
  double SoundSpeed[MeshSize+1];   //sound speed at cell edges.
  double uE[MeshSize+1];            // velocity at cell edges.
  int j;
  for (int i = 0; i <= MeshSize; i++) {
    j = num_ghost + i -1;
    SoundSpeed[i] = 0.5*( sqrt(gamma*P_tomb[j]/rho_tomb[j]) + sqrt(gamma*P_tomb[j+1]/rho_tomb[j+1]) );
    uE[i] = 0.5*(u_tomb[num_ghost+i-1] + u_tomb[num_ghost+i]);
    *(ptr_lambda_max+i) = abs(uE[i]) + SoundSpeed[i];
  }
}
 
void Add_ArtiDissipation(int MeshSize, int num_ghost, double lambda_max[], double P_tomb[], \
			 double Consrv1_tomb[], double Consrv2_tomb[], double Consrv3_tomb[], \
			 double *ptr_Flux1, double *ptr_Flux2, double *ptr_Flux3)
{
  double K2 = 0.25;    // 0.25 <= K2 <= 0.5;
  double K4 = 1.0/64;   // 1/64 <= K4 <= 1/32;
  double dis[MeshSize+1];
  double eps2[MeshSize+1], eps4[MeshSize+1];
  double D1U1[MeshSize+1], D1U2[MeshSize+1], D1U3[MeshSize+1];
  double D3U1[MeshSize+1], D3U2[MeshSize+1], D3U3[MeshSize+1];
  double nu[MeshSize+2*num_ghost];
  int j;
  
  for (int i = 1; i < MeshSize+2*num_ghost-1; i++)
    {
      nu[i] = abs( (P_tomb[i+1] -2*P_tomb[i] + P_tomb[i-1])/(P_tomb[i+1] + 2*P_tomb[i] + P_tomb[i-1]) );
    }
  for (int i = 0; i < MeshSize+1; i++)
    {
      j = num_ghost+i-1;
      eps2[i] = Max(nu[j-1],nu[j]);
      eps2[i] = Max(eps2[i],nu[j+1]);
      eps2[i] = Max(eps2[i],nu[j+2]);
      eps2[i] = eps2[i]*K2;
      
      eps4[i] = Max(0.0, K4-eps2[i]);
    }
  for (int i =0; i < MeshSize+1; i++)
    {
      j = num_ghost + i -1;

      D1U1[i] = lambda_max[i]*eps2[i]*(Consrv1_tomb[j+1]-Consrv1_tomb[j]);
      D1U2[i] = lambda_max[i]*eps2[i]*(Consrv2_tomb[j+1]-Consrv2_tomb[j]);
      D1U3[i] = lambda_max[i]*eps2[i]*(Consrv3_tomb[j+1]-Consrv3_tomb[j]);

      D3U1[i] = lambda_max[i]*eps4[i]*(Consrv1_tomb[j+2]-3*Consrv1_tomb[j+1]+3*Consrv1_tomb[j]-Consrv1_tomb[j-1]);
      D3U2[i] = lambda_max[i]*eps4[i]*(Consrv2_tomb[j+2]-3*Consrv2_tomb[j+1]+3*Consrv2_tomb[j]-Consrv2_tomb[j-1]);
      D3U3[i] = lambda_max[i]*eps4[i]*(Consrv3_tomb[j+2]-3*Consrv3_tomb[j+1]+3*Consrv3_tomb[j]-Consrv3_tomb[j-1]);

      *(ptr_Flux1+i) = *(ptr_Flux1+i) -D1U1[i] + D3U1[i];
      *(ptr_Flux2+i) = *(ptr_Flux2+i) -D1U2[i] + D3U2[i];
      *(ptr_Flux3+i) = *(ptr_Flux3+i) -D1U3[i] + D3U3[i];
    }
}

void Calc_S(int MeshSize, double P[], double dAcdx[], double *ptr_S1, double *ptr_S2, double *ptr_S3)
{
  for (int i=0; i<MeshSize; i++)
    {
      *(ptr_S1+i) = 0.0;
      *(ptr_S2+i) = P[i]*dAcdx[i];
      *(ptr_S3+i) = 0.0;
    }
}

double Max(double a, double b)
{
  return ( 0.5*(a+b + abs(a-b)) );
}

double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}
 
void InflowBC(int MeshSize, int num_ghost, double u_tomb[], double P0, double T0, double gamma, double R, \
	      double *ptr_rho_tomb, double *ptr_P_tomb)
{
  double phi;
  for (int i=0; i<num_ghost; i++)
    {
      phi = Calc_u2phi(u_tomb[i], T0, gamma, R);
      *(ptr_P_tomb+i) = Calc_phi2P(phi, P0, gamma);
      *(ptr_rho_tomb+i) = Calc_phi2rho(phi, P0, T0, R, gamma);
    }
}

void OutflowBC(int flag_isen, int MeshSize, int num_ghost, double Pback, double *ptr_P_tomb)  
{
  switch (flag_isen)
    {
    case 0:              // Normal shock (subsonic outflow)
      for (int i=MeshSize+num_ghost; i<MeshSize+2*num_ghost; i++)
	{
	  *(ptr_P_tomb+i) = 2*Pback - *(ptr_P_tomb+i-1);
	}
      break;

    case 1:              // Isentropic (supersonic outflow)
      break;
    }
}

void Calc_Res(int MeshSize, double Flux[], double Area[], double Source[], double dx[], double *ptr_Res)
{
  for (int i=0; i<MeshSize; i++)
    {
      *(ptr_Res+i) = Flux[i+1]*Area[i+1] - Flux[i]*Area[i] - Source[i]*dx[i];
    }
}

double Calc_L2Norm(int MeshSize, double X[])
{
  double sum = 0;
  for (int i=0; i<MeshSize; i++)
    {
      sum = sum + X[i]*X[i];
    }
  return sqrt(sum/MeshSize);
}

void Calc_upwind_PrimLR(int MeshSize, int num_ghost, double upwind_eps, double upwind_kappa,
			    double prm_tomb[], double *ptr_prm_L, double *ptr_prm_R)
{
  int j;
  double denom;
  double eps = 1e-6;
  double r1_pos, r2_pos, r2_neg, r3_neg;
  double phi1_pos, phi2_pos, phi2_neg, phi3_neg;
  
  for (int i=0; i<MeshSize+1; i++) {
    j = i + num_ghost;

    denom = (prm_tomb[j-1]-prm_tomb[j-2]);
    denom = SIGN( Max( abs(denom), eps ) , denom );
    r1_pos = (prm_tomb[j]-prm_tomb[j-1])/denom;      //  r2[i-1/2]_pos

    denom = (prm_tomb[j]-prm_tomb[j-1]);
    denom = SIGN( Max( abs(denom), eps ) , denom );
    r2_pos = (prm_tomb[j+1]-prm_tomb[j])/denom;       //  r2[i+1/2]_pos
    
    denom = (prm_tomb[j]-prm_tomb[j-1]);
    denom = SIGN( Max( abs(denom), eps ) , denom );
    r2_neg = (prm_tomb[j-1]-prm_tomb[j-2])/denom;     //  r2[i+1/2]_neg

    denom = (prm_tomb[j+1]-prm_tomb[j]);
    denom = SIGN( Max( abs(denom), eps ) , denom );
    r3_neg = (prm_tomb[j]-prm_tomb[j-1])/denom;       //  r3[i+3/2]_neg

    // no Limiters
    /*
    phi1_pos = 1.0;
    phi2_pos = 1.0;
    phi2_neg = 1.0;
    phi3_neg = 1.0;
    */
    
    // using Van Leer Limiters
    /*
    phi1_pos = (r1_pos + abs(r1_pos)) / (1 + r1_pos);
    phi2_pos = (r2_pos + abs(r2_pos)) / (1 + r2_pos);
    phi2_neg = (r2_neg + abs(r2_neg)) / (1 + r2_neg);
    phi3_neg = (r3_neg + abs(r3_neg)) / (1 + r3_neg);
    */

    // using Van Albada Limiters
    
    if (r1_pos < 0.0)
      phi1_pos = 0.0;
    else
      phi1_pos = (r1_pos + r1_pos*r1_pos) / (1 + r1_pos*r1_pos);
    
    if (r2_pos < 0.0)
      phi2_pos = 0.0;
    else
      phi2_pos = (r2_pos + r2_pos*r2_pos) / (1 + r2_pos*r2_pos);

    if (r2_neg < 0.0)
      phi2_neg = 0.0;
    else
    phi2_neg = (r2_neg + r2_neg*r2_neg) / (1 + r2_neg*r2_neg);

    if (r3_neg < 0.0)
      phi3_neg = 0.0;
    else
    phi3_neg = (r3_neg + r3_neg*r3_neg) / (1 + r3_neg*r3_neg);
    
    
    ptr_prm_L[i] = prm_tomb[j-1] + 0.25*upwind_eps*((1-upwind_kappa) * phi1_pos * (prm_tomb[j-1]-prm_tomb[j-2]) \
						    + (1+upwind_kappa) * phi2_neg * (prm_tomb[j]-prm_tomb[j-1]));
    ptr_prm_R[i] = prm_tomb[j] - 0.25*upwind_eps*((1+upwind_kappa) * phi2_pos * (prm_tomb[j]-prm_tomb[j-1]) \
						  + (1-upwind_kappa) * phi3_neg * (prm_tomb[j+1]-prm_tomb[j]));
  }
}

void Calc_upwind_Flux_VL(int MeshSize, double gamma, double R, double rho_L[], double rho_R[], double u_L[], \
			 double u_R[], double P_L[], double P_R[], double *ptr_F1, double *ptr_F2, double *ptr_F3)
{
  double M_L, M_R, M_pos, M_neg;
  double beta_L, beta_R;
  double alpha_pos, alpha_neg;
  double c_pos, c_neg;
  double Vs_L, Vs_R;                      // sound speed
  double D_pos, D_neg;
  double Pbar_pos, Pbar_neg;
  double ht_L, ht_R;
  double Fc1, Fc2, Fc3;
  double Fp1, Fp2, Fp3;
  
  for (int i=0; i<MeshSize+1; i++) {
    Vs_L = sqrt(gamma*P_L[i]/rho_L[i]);
    Vs_R = sqrt(gamma*P_R[i]/rho_R[i]);
    
    M_L = u_L[i]/Vs_L;
    M_pos = 0.25*(M_L + 1)*(M_L + 1);
    M_R = u_R[i]/Vs_R;
    M_neg = -0.25*(M_R - 1)*(M_R - 1);
   
    beta_L = -Max(0, 1 - int(abs(M_L)));
    beta_R = -Max(0, 1 - int(abs(M_R)));

    alpha_pos = 0.5*(1 + SIGN(1,M_L));
    alpha_neg = 0.5*(1 - SIGN(1,M_R));

    c_pos = alpha_pos*(1+beta_L)*M_L - beta_L*M_pos;
    c_neg = alpha_neg*(1+beta_R)*M_R - beta_R*M_neg;

    Pbar_pos = M_pos*(-M_L + 2);
    Pbar_neg = M_neg*(-M_R - 2);

    D_pos = alpha_pos*(1+beta_L) - beta_L*Pbar_pos;
    D_neg = alpha_neg*(1+beta_R) - beta_R*Pbar_neg;

    ht_L = (gamma/(gamma-1))*P_L[i]/rho_L[i]+u_L[i]*u_L[i]/2;
    ht_R = (gamma/(gamma-1))*P_R[i]/rho_R[i]+u_R[i]*u_R[i]/2;
    
    Fc1 = rho_L[i]*Vs_L*c_pos + rho_R[i]*Vs_R*c_neg;
    Fc2 = rho_L[i]*Vs_L*c_pos*u_L[i] + rho_R[i]*Vs_R*c_neg*u_R[i];
    Fc3 = rho_L[i]*Vs_L*c_pos*ht_L + rho_R[i]*Vs_R*c_neg*ht_R;

    Fp1 = 0;
    Fp2 = D_pos*P_L[i] + D_neg*P_R[i];
    Fp3 = 0;

    ptr_F1[i] = Fc1 + Fp1;
    ptr_F2[i] = Fc2 + Fp2;
    ptr_F3[i] = Fc3 + Fp3;
  }  
}


void Calc_upwind_Flux_Roe(int MeshSize, double gamma, double R, double rho_L[], double rho_R[], double u_L[], \
			  double u_R[], double P_L[], double P_R[], double *ptr_F1, double *ptr_F2, double *ptr_F3)
{
  int ne = 3;           // Number of equations we have in our Euler system

  double eps = 0.1;     // the eps for the modified Roe averaged eigenvalues
  double rho_Roe;
  double u_Roe;
  double ht_L, ht_R;
  double ht_Roe;
  double Vs_Roe;        // Roe averaged sound speed
  double R_Roe;         // sqrt(rho_R/rho_L)
  double drho, du, dP;

  double abs_lambda1_Roe, abs_lambda2_Roe, abs_lambda3_Roe;
  double r1_Roe[ne], r2_Roe[ne], r3_Roe[ne];
  double dw1, dw2, dw3;

  double F1_L, F1_R, F2_L, F2_R, F3_L, F3_R;

  for (int i=0; i<MeshSize+1; i++) {
    R_Roe = sqrt(rho_R[i]/rho_L[i]);
    rho_Roe = R_Roe*rho_L[i];
    u_Roe = (R_Roe*u_R[i] + u_L[i]) / (R_Roe + 1);
    ht_L = (gamma/(gamma-1))*P_L[i]/rho_L[i]+u_L[i]*u_L[i]/2;
    ht_R = (gamma/(gamma-1))*P_R[i]/rho_R[i]+u_R[i]*u_R[i]/2;
    ht_Roe = (R_Roe*ht_R + ht_L) / (R_Roe + 1);
    Vs_Roe = sqrt( (gamma-1)*(ht_Roe - u_Roe*u_Roe/2) );
    
    abs_lambda1_Roe = abs(u_Roe);
    abs_lambda2_Roe = abs(u_Roe + Vs_Roe);
    abs_lambda3_Roe = abs(u_Roe - Vs_Roe);
    // calc modified lambda
    if (abs_lambda1_Roe < 2*eps*Vs_Roe)
      abs_lambda1_Roe = abs_lambda1_Roe*abs_lambda1_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;
    if (abs_lambda2_Roe < 2*eps*Vs_Roe)
      abs_lambda2_Roe = abs_lambda2_Roe*abs_lambda2_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;
    if (abs_lambda3_Roe < 2*eps*Vs_Roe)
      abs_lambda3_Roe = abs_lambda3_Roe*abs_lambda3_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;

    r1_Roe[0] = 1.0;
    r1_Roe[1] = u_Roe;
    r1_Roe[2] = u_Roe*u_Roe/2;

    r2_Roe[0] = (rho_Roe/(2*Vs_Roe));
    r2_Roe[1] = (rho_Roe/(2*Vs_Roe))*(u_Roe+Vs_Roe);
    r2_Roe[2] = (rho_Roe/(2*Vs_Roe))*(ht_Roe+u_Roe*Vs_Roe);

    r3_Roe[0] = -(rho_Roe/(2*Vs_Roe));
    r3_Roe[1] = -(rho_Roe/(2*Vs_Roe))*(u_Roe-Vs_Roe);
    r3_Roe[2] = -(rho_Roe/(2*Vs_Roe))*(ht_Roe-u_Roe*Vs_Roe);

    drho = rho_R[i] - rho_L[i];
    du = u_R[i] - u_L[i];
    dP = P_R[i] - P_L[i];

    dw1 = drho - dP/(Vs_Roe*Vs_Roe);
    dw2 = du + dP/(rho_Roe*Vs_Roe);
    dw3 = du - dP/(rho_Roe*Vs_Roe);

    F1_L = rho_L[i]*u_L[i];
    F1_R = rho_R[i]*u_R[i];

    F2_L = rho_L[i]*u_L[i]*u_L[i] + P_L[i];
    F2_R = rho_R[i]*u_R[i]*u_R[i] + P_R[i];

    F3_L = rho_L[i]*u_L[i]*ht_L;
    F3_R = rho_R[i]*u_R[i]*ht_R;

    ptr_F1[i] = 0.5*(F1_L + F1_R) - 0.5*( abs(abs_lambda1_Roe)*dw1*r1_Roe[0] + abs(abs_lambda2_Roe)*dw2*r2_Roe[0] + abs(abs_lambda3_Roe)*dw3*r3_Roe[0]);
    ptr_F2[i] = 0.5*(F2_L + F2_R) - 0.5*( abs(abs_lambda1_Roe)*dw1*r1_Roe[1] + abs(abs_lambda2_Roe)*dw2*r2_Roe[1]
					  + abs(abs_lambda3_Roe)*dw3*r3_Roe[1]);
    ptr_F3[i] = 0.5*(F3_L + F3_R) - 0.5*( abs(abs_lambda1_Roe)*dw1*r1_Roe[2] + abs(abs_lambda2_Roe)*dw2*r2_Roe[2]
					  + abs(abs_lambda3_Roe)*dw3*r3_Roe[2]);     
  }
}

//===============================================================================================================//
//                                                 I/O Functions                                                 //
//===============================================================================================================//
void OutputSingleArray(string FileName, int ArrSize, double x[], double Arr[], int iteration)
{
  ostringstream StrConvert;
  StrConvert << iteration;
  string num_iter = StrConvert.str();
  
  string Address = "./OutputFiles/";
  string Suffix = ".txt";
  Address = Address + FileName + num_iter + Suffix;
  
  ofstream outfile;
  outfile.open(Address);

  for (int i = 0; i <ArrSize; i++)
    {
      outfile << setprecision(14)<< i <<"     " << x[i] <<"     " << Arr[i] << endl;
    }
  outfile.close();
}

void OutputSingleVector(string FileName, int VecSize, vector<double> Vec)
{
  string Address = "./OutputFiles/";
  string Suffix = ".txt";
  Address = Address + FileName + Suffix;
  
  ofstream outfile;
  outfile.open(Address);

  for (int i = 0; i <VecSize; i++)
    {
      outfile << setprecision(14)<< i <<"     "<< Vec[i] << endl;
    }
  outfile.close();
}

void CoutTitle(string TitleName, string LineStyle, string LengthStyle)
{
  int Num;
  int TitleLength = TitleName.size();
  if (LengthStyle == "long")
    Num = 100;
  else if (LengthStyle == "short")
    Num = 50;
  else {
    cout<<"ERROR: Wrong Input for 'LengthStyle' in void function OutputTitle(string TitleName, ";
    cout<<"string LineStyle, string LengthStyle), Please input: 'long' or 'short' \n";
  }
  for (int i=0; i<Num; i++)
    cout<<LineStyle;
  cout<<"\n";
  for (int i=0; i<Num-TitleLength-8; i++) {
    cout<<" ";
    if(i == int( (Num-TitleLength-8) /2) )
      cout<<"*   "<<TitleName<<"   *";
  }
  cout<<"\n";
  for (int i=0; i<Num; i++)
    cout<<LineStyle;
  cout<<"\n";
}

//===============================================================================================================//
//                                                  *   The End   *                                              //
//===============================================================================================================//
