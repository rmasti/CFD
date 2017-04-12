//===============================================================================================================//
//                                        *  AOE6145 CFD Final project  *                                        //
//===============================================================================================================//
// 2D Euler equations solver
// Upwind Schemes ( 1st&2nd order of accuracy, MUSCL extrapolation, Van Leer & Roe, Flux(slope) limiters )
// Method of Manufactured Solution (code verification)
/* Cases: 
   1. Curvilinear structured meshes;
   2. 2D 30 degree inlet;  
   3. NACA64A series airfoil with a thickness-to-chord ratio of 6%
*/

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

int flag_case;                // 1 for curve Linear; 2 for 30 degree inlet; 3 for NACA Airfoil
int flag_mesh;                // chose from 1-4 (1-6 for curve linear case) for meshes from coarse to fine
int flag_upwind;              // 1 for Van Leer; 2 for Roe
int flag_limiter;             // set 0 to turn off limiter, set 1 to turn on limiter
int flag_supersonic;          // 0 for subsonic; 1 for supersonic ( only valid for flag_case == 1 )
int flag_AOA;                 // 1 for 0 deg of AOA, 2 for 8 deg of AOA ( only valid for falg_case == 3 )
int RK_order;                 // order of numerical order for RK_order
double upwind_eps;            
string debug = "true";
string Output_Folder;

#define Pi 3.1415926
#define gamma 1.4
#define R 287.0

//===============================================================================================================//
//                                            Declare Functions                                                  //
//===============================================================================================================//
void Calc_Area(int ni, int nj, vector<vector<double> > x, vector<vector<double> > y, vector<vector<double> > &Ai,
	       vector<vector<double> > &Aj);
void Calc_NormVec(vector<vector<double> > x, vector<vector<double> > y, vector<vector<double> > Ai,
		  vector<vector<double> > Aj, vector<vector<double> > &norm_i_x, vector<vector<double> > &norm_i_y,
		  vector<vector<double> > &norm_j_x, vector<vector<double> > &norm_j_y);
void Calc_Volumn(int ni, int nj, vector<vector<double> > x, vector<vector<double> > y,
		 vector<vector<double> > &Volumn);
void Resize_2dVec(int size_i, int size_j, vector<vector<double> > &Vec2d);
double Vector_Cross_norm(double x1, double y1, double x2, double y2);
void GoTomb_2d_extrapolate(int num_ghost, vector<vector<double> > Vec2d, vector<vector<double> > &Vec2d_tomb);
void GoTomb_2d(int num_ghost, vector<vector<double> > Vec2d, vector<vector<double> > &Vec2d_tomb);
void OutTomb_2d(int num_ghost, vector<vector<double> > Vec2d_tomb, vector<vector<double> > &Vec2d);
void Ini_rho(vector<vector<double> > xc, vector<vector<double> > yc, vector<vector<double> > &rho);
void Ini_u(vector<vector<double> > xc, vector<vector<double> > yc,
	   vector<vector<double> > &u);
void Ini_v(vector<vector<double> > xc, vector<vector<double> > yc,
	   vector<vector<double> > &v);
void Ini_P(vector<vector<double> > xc, vector<vector<double> > yc, vector<vector<double> > &P);
void Calc_Temperature(vector<vector<double> > rho, vector<vector<double> > P, vector<vector<double> > &T);
void solve_MMS(vector<vector<double> > xc, vector<vector<double> > yc, vector<vector<double> > &rho,
	       vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &P);
void solve_MMS_source(vector<vector<double> > xc, vector<vector<double> > yc,
		      vector<vector<double> > &S1, vector<vector<double> > &S2, vector<vector<double> > &S3,
		      vector<vector<double> > &S4);
void Enforce_BC_MMS(int num_ghost, vector<vector<double> > prm_tomb_MMS, vector<vector<double> > &prm_tomb);
void Enforce_BC(int num_ghost, vector<vector<double> > ni_x, vector<vector<double> > ni_y,
		vector<vector<double> > nj_x, vector<vector<double> > nj_y,  vector<vector<double> > &rho_tomb,
		vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb,
		vector<vector<double> > &P_tomb, vector<vector<double> > &T_tomb);
void Symmetric_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
		   vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb,
		   vector<vector<double> > &P_tomb);
void Inlet_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
	      vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb);
void Outlet_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
	       vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb);
void SlipWall_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > ni_x,
		 vector<vector<double> > ni_y, vector<vector<double> > nj_x, vector<vector<double> > nj_y,
		 vector<vector<double> > &rho_tomb, vector<vector<double> > &u_tomb,
		 vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb, vector<vector<double> > &T_tomb);
void FarField_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
		 vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb);
void Periodic_BC(int A_Begin[], int A_End[], int B_Begin[], int B_End[], int num_ghost,
		 vector<vector<double> > &rho_tomb, vector<vector<double> > &u_tomb,
		 vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb);
void Switch_Prm2Consrv(vector<vector<double> > V1, vector<vector<double> > V2,
		       vector<vector<double> > V3, vector<vector<double> > V4, vector<vector<double> > &U1,
		       vector<vector<double> > &U2, vector<vector<double> > &U3, vector<vector<double> > &U4);
void Switch_Consrv2Prm(vector<vector<double> > U1, vector<vector<double> > U2,
		       vector<vector<double> > U3, vector<vector<double> > U4, vector<vector<double> > &V1,
		       vector<vector<double> > &V2, vector<vector<double> > &V3, vector<vector<double> > &V4);
void Calc_upwind_Prim_i(int num_ghost, double upwind_eps, double upwind_kappa, vector<vector<double> > prm_tomb,
			vector<vector<double> > &prm_Left, vector<vector<double> > &prm_Right);
void Calc_upwind_Prim_j(int num_ghost, double upwind_eps, double upwind_kappa, vector<vector<double> > prm_tomb,
			vector<vector<double> > &prm_Lower, vector<vector<double> > &prm_Upper);
void MUSCL_prm(int num_ghost, double upwind_eps, double upwind_kappa, vector<vector<double> > rho_tomb,
	       vector<vector<double> > u_tomb, vector<vector<double> > v_tomb, vector<vector<double> > P_tomb,
	       vector<vector<double> > &rho_Left, vector<vector<double> > &rho_Right,
	       vector<vector<double> > &rho_Lower, vector<vector<double> > &rho_Upper,
	       vector<vector<double> > &u_Left, vector<vector<double> > &u_Right,
	       vector<vector<double> > &u_Lower, vector<vector<double> > &u_Upper,
	       vector<vector<double> > &v_Left, vector<vector<double> > &v_Right,
	       vector<vector<double> > &v_Lower, vector<vector<double> > &v_Upper,
	       vector<vector<double> > &P_Left, vector<vector<double> > &P_Right,
	       vector<vector<double> > &P_Lower, vector<vector<double> > &P_Upper);
void Calc_Flux_VL(vector<vector<double> > nx, vector<vector<double> > ny,
		  vector<vector<double> > rho_L, vector<vector<double> > rho_R, vector<vector<double> > u_L,
		  vector<vector<double> > u_R, vector<vector<double> > v_L, vector<vector<double> > v_R,
		  vector<vector<double> > P_L, vector<vector<double> > P_R, vector<vector<double> > &F1,
		  vector<vector<double> > &F2, vector<vector<double> > &F3, vector<vector<double> > &F4);
void Calc_Flux_Roe(vector<vector<double> > nx, vector<vector<double> > ny,
		   vector<vector<double> > rho_L, vector<vector<double> > rho_R, vector<vector<double> > u_L,
		   vector<vector<double> > u_R, vector<vector<double> > v_L, vector<vector<double> > v_R,
		   vector<vector<double> > P_L, vector<vector<double> > P_R, vector<vector<double> > &F1,
		   vector<vector<double> > &F2, vector<vector<double> > &F3, vector<vector<double> > &F4);
void Calc_Flux(vector<vector<double> > nx, vector<vector<double> > ny,
	       vector<vector<double> > rho_L, vector<vector<double> > rho_R, vector<vector<double> > u_L,
	       vector<vector<double> > u_R, vector<vector<double> > v_L, vector<vector<double> > v_R,
	       vector<vector<double> > P_L, vector<vector<double> > P_R, vector<vector<double> > &F1,
	       vector<vector<double> > &F2, vector<vector<double> > &F3, vector<vector<double> > &F4);
void Calc_Vs(vector<vector<double> > rho, vector<vector<double> > P, vector<vector<double> > &Vs);
void Calc_dt(double CFL, vector<vector<double> > Volumn, vector<vector<double> > Ai,
	     vector<vector<double> > Aj, vector<vector<double> > ni_x, vector<vector<double> > ni_y,
	     vector<vector<double> > nj_x, vector<vector<double> > nj_y, vector<vector<double> > u,
	     vector<vector<double> > v, vector<vector<double> > Vs, vector<vector<double> > &dt);
void Calc_Res(vector<vector<double> > Fi, vector<vector<double> > Fj, vector<vector<double> > S,
	      vector<vector<double> > Ai, vector<vector<double> > Aj, vector<vector<double> > volumn,
	      vector<vector<double> > &Res);
void Calc_Error(vector<vector<double> > Simulation, vector<vector<double> > Exact, vector<vector<double> > &Error);
void RK_marching(int stage, vector<vector<double> > Consrv_old, double dt_global, vector<vector<double> > Volumn,
		 vector<vector<double> > Res, vector<vector<double> > &Consrv);
double Calc_L2Norm(vector<vector<double> > X);
double Max(double a, double b);
double Min(double a, double b);
double Min_Element(vector<vector<double> > Vec);
double SIGN(double a, double b);

//--------------------------------------------- * I/O Functions * -----------------------------------------------//
void CoutTitle(string TitleName, string LineStyle, string LengthStyle);
double load_inputNum(string FileName, string Var);
string read_mesh_name();
void InputMesh(int &ni, int &nj, int &nk, vector<vector<double> > &x, vector<vector<double> > &xc,
	       vector<vector<double> > &y, vector<vector<double> > &yc, vector<vector<double> > &z,
	       vector<vector<double> > &zc);
void build_case_folder();
void Output_2D_Vector(string FileName, vector<vector<double> > Vec2d, int iteration);
void Output_1D_Vector(string FileName, vector<double> Vec);
void Save_data_to_folder();

//===============================================================================================================//
//                                              Main Functions                                                   //
//===============================================================================================================//
int main(int argc, char *argv[])
{
  //***********************************************************************************************************//
  //                                         Input Variables and Mesh                                          //
  //***********************************************************************************************************//
  const char * InputFile = argv[1];
  flag_case = load_inputNum(InputFile, "flag_case");
  flag_mesh = load_inputNum(InputFile, "flag_mesh");
  flag_supersonic = load_inputNum(InputFile, "flag_supersonic");
  flag_AOA = load_inputNum(InputFile, "flag_AOA");
  flag_upwind = load_inputNum(InputFile, "flag_upwind");
  flag_limiter = load_inputNum(InputFile, "flag_limiter");
  RK_order = load_inputNum(InputFile, "RK_order");
  const double CFL = load_inputNum(InputFile, "CFL");
  const int max_iter = load_inputNum(InputFile, "max_iter");
  const double Final_res = load_inputNum(InputFile, "Final_res");
  const int WriteInterval = load_inputNum(InputFile, "WriteInterval");
  const int PrintInterval = load_inputNum(InputFile, "PrintInterval");
  upwind_eps = load_inputNum(InputFile, "upwind_eps");
  double upwind_kappa = load_inputNum(InputFile, "upwind_kappa");

  build_case_folder();
  
  const int num_ghost = 3;                   // number of ghost cells on each side;

  int iteration = 0;
  
  int ni;   // numer of cells on the i direction
  int nj;   // numer of cells on the j direction
  int nk;   // numer of cells on the k direction
  vector<vector<double> > x;
  vector<vector<double> > y;
  vector<vector<double> > z;
  vector<vector<double> > xc;
  vector<vector<double> > yc;
  vector<vector<double> > zc;
  
  InputMesh(ni,nj,nk,x,xc,y,yc,z,zc);
  
  CoutTitle("Input Variables", "=", "short");
  cout<<" * flag_case:          "<<flag_case<<endl;
  cout<<" * flag_mesh:          "<<flag_mesh<<endl;
  if (flag_case == 1)
    cout<<" * flag_supersonic:    "<<flag_supersonic<<endl;
  if (flag_case == 3)
    cout<<" * flag_AOA:           "<<flag_AOA<<endl;
  cout<<" * flag_upwind:        "<<flag_upwind<<endl;
  cout<<" * flag_limiter:       "<<flag_limiter<<endl;
  cout<<" * RK_order:           "<<RK_order<<endl;
  cout<<" * CFL:                "<<CFL<<endl;
  cout<<" * max # of iteration: "<<max_iter<<endl;
  cout<<" * final residual:     "<<Final_res<<endl;
  cout<<" * Write Interval:     "<<WriteInterval<<endl;
  cout<<" * Upwind epsilon:     "<<upwind_eps<<endl;
  cout<<" * Upwind kappa:       "<<upwind_kappa<<endl;

  //***********************************************************************************************************//
  //                                         Setting Other Variables                                           //
  //***********************************************************************************************************//
  vector<vector<double> > Ai(ni+1, vector<double>(nj));
  vector<vector<double> > Aj(ni, vector<double>(nj+1));
  vector<vector<double> > norm_i_x(ni+1, vector<double>(nj));
  vector<vector<double> > norm_i_y(ni+1, vector<double>(nj));
  vector<vector<double> > norm_j_x(ni, vector<double>(nj+1));
  vector<vector<double> > norm_j_y(ni, vector<double>(nj+1));
  
  vector<vector<double> > Volumn(ni, vector<double>(nj));
  
  vector<vector<double> > rho(ni, vector<double>(nj));
  vector<vector<double> > u(ni, vector<double>(nj));
  vector<vector<double> > v(ni, vector<double>(nj));
  vector<vector<double> > P(ni, vector<double>(nj));

  vector<vector<double> > rho_MMS(ni, vector<double>(nj));
  vector<vector<double> > u_MMS(ni, vector<double>(nj));
  vector<vector<double> > v_MMS(ni, vector<double>(nj));
  vector<vector<double> > P_MMS(ni, vector<double>(nj));

  vector<vector<double> > Vs(ni, vector<double>(nj));                    // speed of sound

  vector<vector<double> > U1(ni, vector<double>(nj));
  vector<vector<double> > U2(ni, vector<double>(nj));
  vector<vector<double> > U3(ni, vector<double>(nj));
  vector<vector<double> > U4(ni, vector<double>(nj));
  
  vector<vector<double> > U1_RK(ni, vector<double>(nj));
  vector<vector<double> > U2_RK(ni, vector<double>(nj));
  vector<vector<double> > U3_RK(ni, vector<double>(nj));
  vector<vector<double> > U4_RK(ni, vector<double>(nj));

  vector<vector<double> > Fi_1(ni+1, vector<double>(nj));
  vector<vector<double> > Fi_2(ni+1, vector<double>(nj));
  vector<vector<double> > Fi_3(ni+1, vector<double>(nj));
  vector<vector<double> > Fi_4(ni+1, vector<double>(nj));
  
  vector<vector<double> > Fj_1(ni, vector<double>(nj+1));
  vector<vector<double> > Fj_2(ni, vector<double>(nj+1));
  vector<vector<double> > Fj_3(ni, vector<double>(nj+1));
  vector<vector<double> > Fj_4(ni, vector<double>(nj+1));

  vector<vector<double> > S1(ni, vector<double>(nj));
  vector<vector<double> > S2(ni, vector<double>(nj));
  vector<vector<double> > S3(ni, vector<double>(nj));
  vector<vector<double> > S4(ni, vector<double>(nj));

  vector<vector<double> > Res1(ni, vector<double>(nj));
  vector<vector<double> > Res2(ni, vector<double>(nj));
  vector<vector<double> > Res3(ni, vector<double>(nj));
  vector<vector<double> > Res4(ni, vector<double>(nj));

  vector<vector<double> > Error1(ni, vector<double>(nj));
  vector<vector<double> > Error2(ni, vector<double>(nj));
  vector<vector<double> > Error3(ni, vector<double>(nj));
  vector<vector<double> > Error4(ni, vector<double>(nj));

  vector<vector<double> > dt(ni, vector<double>(nj));
  double dt_global;

  double L2norm1_initial;
  double L2norm2_initial;
  double L2norm3_initial;
  double L2norm4_initial;

  double L2norm1;
  double L2norm2;
  double L2norm3;
  double L2norm4;

  double iterE_rho_initial;
  double iterE_u_initial;
  double iterE_v_initial;
  double iterE_P_initial;
  
  double iterE_rho;
  double iterE_u;
  double iterE_v;
  double iterE_P;

  vector<double> L2norm1_rel;
  vector<double> L2norm2_rel;
  vector<double> L2norm3_rel;
  vector<double> L2norm4_rel;

  vector<double> iterError_rho;
  vector<double> iterError_u;
  vector<double> iterError_v;
  vector<double> iterError_P;

  //***********************************************************************************************************//
  //                                      Defined Upwind Scheme Variables                                      //
  //***********************************************************************************************************//
  vector<vector<double> > rho_Left(ni+1, vector<double>(nj));
  vector<vector<double> > u_Left(ni+1, vector<double>(nj));
  vector<vector<double> > v_Left(ni+1, vector<double>(nj));
  vector<vector<double> > P_Left(ni+1, vector<double>(nj));

  vector<vector<double> > rho_Right(ni+1, vector<double>(nj));
  vector<vector<double> > u_Right(ni+1, vector<double>(nj));
  vector<vector<double> > v_Right(ni+1, vector<double>(nj));
  vector<vector<double> > P_Right(ni+1, vector<double>(nj));

  vector<vector<double> > rho_Lower(ni, vector<double>(nj+1));
  vector<vector<double> > u_Lower(ni, vector<double>(nj+1));
  vector<vector<double> > v_Lower(ni, vector<double>(nj+1));
  vector<vector<double> > P_Lower(ni, vector<double>(nj+1));
  
  vector<vector<double> > rho_Upper(ni, vector<double>(nj+1));
  vector<vector<double> > u_Upper(ni, vector<double>(nj+1));
  vector<vector<double> > v_Upper(ni, vector<double>(nj+1));
  vector<vector<double> > P_Upper(ni, vector<double>(nj+1));
  
  //***********************************************************************************************************//
  //              Define ghost and tomb Variables (tomb contains both the ACTUAL and GHOST cells)              //
  //***********************************************************************************************************//
  /*
  vector<vector<double> > xc_ghost_i(2*num_ghost, vector<double>(2*num_ghost));
  vector<vector<double> > xc_ghost_j(2*num_ghost, vector<double>(2*num_ghost));
  vector<vector<double> > yc_ghost_i(2*num_ghost, vector<double>(2*num_ghost));
  vector<vector<double> > yc_ghost_j(2*num_ghost, vector<double>(2*num_ghost));
  */
  
  vector<vector<double> > xc_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > yc_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  
  vector<vector<double> > rho_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > u_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > v_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > P_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > T_tomb(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  
  vector<vector<double> > rho_tomb_MMS(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > u_tomb_MMS(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > v_tomb_MMS(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  vector<vector<double> > P_tomb_MMS(ni+2*num_ghost, vector<double>(nj+2*num_ghost));
  
  //***********************************************************************************************************//
  //                                              Initialization                                               //
  //***********************************************************************************************************//
  Calc_Area(ni+1, nj+1, x, y, Ai, Aj);
  Calc_NormVec(x, y, Ai, Aj, norm_i_x, norm_i_y, norm_j_x, norm_j_y);
  Calc_Volumn(ni+1, nj+1, x, y, Volumn);
  GoTomb_2d_extrapolate(num_ghost, xc, xc_tomb);                
  GoTomb_2d_extrapolate(num_ghost, yc, yc_tomb);   

  Ini_rho(xc_tomb, yc_tomb, rho_tomb);
  Ini_u(xc_tomb, yc_tomb, u_tomb);
  Ini_v(xc_tomb, yc_tomb, v_tomb);
  Ini_P(xc_tomb, yc_tomb, P_tomb);
  Calc_Temperature(rho_tomb, P_tomb, T_tomb);

  if (flag_case == 1) {  // if case 1
    solve_MMS_source(xc, yc, S1, S2, S3, S4);                 // calculate MMS source term
    solve_MMS(xc_tomb, yc_tomb, rho_tomb_MMS, u_tomb_MMS, v_tomb_MMS, P_tomb_MMS);   // serve as exact solution
    solve_MMS(xc, yc, rho_MMS, u_MMS, v_MMS, P_MMS);
    
    Enforce_BC_MMS(num_ghost, rho_tomb_MMS, rho_tomb);
    Enforce_BC_MMS(num_ghost, u_tomb_MMS, u_tomb);
    Enforce_BC_MMS(num_ghost, v_tomb_MMS, v_tomb);
    Enforce_BC_MMS(num_ghost, P_tomb_MMS, P_tomb);
  }
  
  else {  // if case 2 or case 3
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nj; j++) {
	S1[i][j] = 0.0;
	S2[i][j] = 0.0;
	S3[i][j] = 0.0;
	S4[i][j] = 0.0;
      }
    }
    Enforce_BC(num_ghost, norm_i_x, norm_i_y, norm_j_x, norm_j_y, rho_tomb, u_tomb, v_tomb, P_tomb, T_tomb);
  }
 
  MUSCL_prm(num_ghost, upwind_eps, upwind_kappa, rho_tomb, u_tomb, v_tomb, P_tomb, rho_Left, rho_Right,
	    rho_Lower, rho_Upper, u_Left, u_Right, u_Lower, u_Upper, v_Left, v_Right, v_Lower, v_Upper,
	    P_Left, P_Right, P_Lower, P_Upper);
  
  Calc_Flux(norm_i_x, norm_i_y, rho_Left, rho_Right, u_Left, u_Right, v_Left, v_Right,
	    P_Left, P_Right, Fi_1, Fi_2, Fi_3, Fi_4);
  Calc_Flux(norm_j_x, norm_j_y, rho_Lower, rho_Upper, u_Lower, u_Upper, v_Lower, v_Upper,
	    P_Lower, P_Upper, Fj_1, Fj_2, Fj_3, Fj_4);

  OutTomb_2d(num_ghost, rho_tomb, rho);
  OutTomb_2d(num_ghost, u_tomb, u);
  OutTomb_2d(num_ghost, v_tomb, v);
  OutTomb_2d(num_ghost, P_tomb, P);

  Switch_Prm2Consrv(rho, u, v, P, U1, U2, U3, U4);
  
  Calc_Res(Fi_1, Fj_1, S1, Ai, Aj, Volumn, Res1);
  Calc_Res(Fi_2, Fj_2, S2, Ai, Aj, Volumn, Res2);
  Calc_Res(Fi_3, Fj_3, S3, Ai, Aj, Volumn, Res3);
  Calc_Res(Fi_4, Fj_4, S4, Ai, Aj, Volumn, Res4);

  L2norm1_initial = Calc_L2Norm(Res1);
  L2norm2_initial = Calc_L2Norm(Res2);
  L2norm3_initial = Calc_L2Norm(Res3);
  L2norm4_initial = Calc_L2Norm(Res4);

  L2norm1_rel.push_back(L2norm1_initial/L2norm1_initial);
  L2norm2_rel.push_back(L2norm2_initial/L2norm2_initial);
  L2norm3_rel.push_back(L2norm3_initial/L2norm3_initial);
  L2norm4_rel.push_back(L2norm4_initial/L2norm4_initial);
  
  if (flag_case ==1) {
    Calc_Error(rho, rho_MMS, Error1);
    Calc_Error(u, u_MMS, Error2);
    Calc_Error(v, v_MMS, Error3);
    Calc_Error(P, P_MMS, Error4);
    
    iterE_rho_initial = Calc_L2Norm(Error1);
    iterE_u_initial = Calc_L2Norm(Error2);
    iterE_v_initial = Calc_L2Norm(Error3);
    iterE_P_initial = Calc_L2Norm(Error4);

    /*
    iterError_rho.push_back(iterE_rho_initial/iterE_rho_initial);
    iterError_u.push_back(iterE_u_initial/iterE_u_initial);
    iterError_v.push_back(iterE_v_initial/iterE_v_initial);
    iterError_P.push_back(iterE_P_initial/iterE_P_initial);
    */

    iterError_rho.push_back(iterE_rho_initial);
    iterError_u.push_back(iterE_u_initial);
    iterError_v.push_back(iterE_v_initial);
    iterError_P.push_back(iterE_P_initial);
  }
  
  Calc_Vs(rho, P, Vs);
  Calc_dt(CFL, Volumn, Ai, Aj, norm_i_x, norm_i_y, norm_j_x, norm_j_y, u, v, Vs, dt);
    
  //-------------------------------------- * Output Initial Variables * ---------------------------------------//
  Output_2D_Vector("x", x, 0);
  Output_2D_Vector("y", y, 0);
  Output_2D_Vector("z", z, 0);
  Output_2D_Vector("xc", xc, 0);
  Output_2D_Vector("yc", yc, 0);
  Output_2D_Vector("zc", zc, 0);
  Output_2D_Vector("Area_i", Ai, 0);
  Output_2D_Vector("Area_j", Aj, 0);
  Output_2D_Vector("norm_i_x", norm_i_x, 0);
  Output_2D_Vector("norm_i_y", norm_i_y, 0);
  Output_2D_Vector("norm_j_x", norm_j_x, 0);
  Output_2D_Vector("norm_j_y", norm_j_y, 0);
  Output_2D_Vector("Volumn", Volumn, 0);
  Output_2D_Vector("xc_tomb", xc_tomb, 0);
  Output_2D_Vector("yc_tomb", yc_tomb, 0);
  
  Output_2D_Vector("rho_tomb", rho_tomb, 0);
  Output_2D_Vector("u_tomb", u_tomb, 0);
  Output_2D_Vector("v_tomb", v_tomb, 0);
  Output_2D_Vector("P_tomb", P_tomb, 0);

  Output_2D_Vector("rho", rho, iteration);
  Output_2D_Vector("u", u, iteration);
  Output_2D_Vector("v", v, iteration);
  Output_2D_Vector("P", P, iteration); 

  Output_2D_Vector("U1", U1, 0);
  Output_2D_Vector("U2", U2, 0);
  Output_2D_Vector("U3", U3, 0);
  Output_2D_Vector("U4", U4, 0);

  Output_2D_Vector("Fi_1", Fi_1, 0);
  Output_2D_Vector("Fi_2", Fi_2, 0);
  Output_2D_Vector("Fi_3", Fi_3, 0);
  Output_2D_Vector("Fi_4", Fi_4, 0);
  Output_2D_Vector("Fj_1", Fj_1, 0);
  Output_2D_Vector("Fj_2", Fj_2, 0);
  Output_2D_Vector("Fj_3", Fj_3, 0);
  Output_2D_Vector("Fj_4", Fj_4, 0);

  Output_2D_Vector("Res1", Res1, 0);
  Output_2D_Vector("Res2", Res2, 0);
  Output_2D_Vector("Res3", Res3, 0);
  Output_2D_Vector("Res4", Res4, 0);
  
  Output_2D_Vector("dt", dt, iteration);

  if (flag_case ==1) {
    Output_2D_Vector("rho_MMS", rho_MMS, 0);
    Output_2D_Vector("u_MMS", u_MMS, 0);
    Output_2D_Vector("v_MMS", v_MMS, 0);
    Output_2D_Vector("P_MMS", P_MMS, 0);
    Output_2D_Vector("S1", S1, 0);
    Output_2D_Vector("S2", S2, 0);
    Output_2D_Vector("S3", S3, 0);
    Output_2D_Vector("S4", S4, 0);
    Output_2D_Vector("Error1", Error1, 0);
    Output_2D_Vector("Error2", Error2, 0);
    Output_2D_Vector("Error3", Error3, 0);
    Output_2D_Vector("Error4", Error4, 0);
  }
  
  //===== Debug =====//
  if (true) {
    Output_2D_Vector("rho_Left", rho_Left, 0);
    Output_2D_Vector("u_Left", u_Left, 0);
    Output_2D_Vector("v_Left", v_Left, 0);
    Output_2D_Vector("P_Left", P_Left, 0);
    Output_2D_Vector("rho_Right", rho_Right, 0);
    Output_2D_Vector("u_Right", u_Right, 0);
    Output_2D_Vector("v_Right", v_Right, 0);
    Output_2D_Vector("P_Right", P_Right, 0);
    Output_2D_Vector("rho_Lower", rho_Lower, 0);
    Output_2D_Vector("u_Lower", u_Lower, 0);
    Output_2D_Vector("v_Lower", v_Lower, 0);
    Output_2D_Vector("P_Lower", P_Lower, 0);
    Output_2D_Vector("rho_Upper", rho_Upper, 0);
    Output_2D_Vector("u_Upper", u_Upper, 0);
    Output_2D_Vector("v_Upper", v_Upper, 0);
    Output_2D_Vector("P_Upper", P_Upper, 0);
  }
  
  //***********************************************************************************************************//
  //                                              Time Marching                                                //
  //***********************************************************************************************************//
  while (iteration < max_iter) {
    iteration++;
    dt_global = Min_Element(dt);
    //dt_global = 1.0e-6;
    //------------------------------- * Runge-Kutta time marching START * -----------------------------------//
    for (int k=0; k<RK_order; k++) {
      RK_marching(k, U1, dt_global, Volumn, Res1, U1_RK);
      RK_marching(k, U2, dt_global, Volumn, Res2, U2_RK);
      RK_marching(k, U3, dt_global, Volumn, Res3, U3_RK);
      RK_marching(k, U4, dt_global, Volumn, Res4, U4_RK);

      Switch_Consrv2Prm(U1_RK, U2_RK, U3_RK, U4_RK, rho, u, v, P);

      GoTomb_2d(num_ghost, rho, rho_tomb);
      GoTomb_2d(num_ghost, u, u_tomb);
      GoTomb_2d(num_ghost, v, v_tomb);
      GoTomb_2d(num_ghost, P, P_tomb);
      Calc_Temperature(rho_tomb, P_tomb, T_tomb);

      if (flag_case != 1) {
	Enforce_BC(num_ghost, norm_i_x, norm_i_y, norm_j_x, norm_j_y, rho_tomb, u_tomb, v_tomb, P_tomb, T_tomb);   
      }
         
      MUSCL_prm(num_ghost, upwind_eps, upwind_kappa, rho_tomb, u_tomb, v_tomb, P_tomb, rho_Left, rho_Right,
		rho_Lower, rho_Upper, u_Left, u_Right, u_Lower, u_Upper, v_Left, v_Right, v_Lower, v_Upper,
		P_Left, P_Right, P_Lower, P_Upper);
      
      Calc_Flux(norm_i_x, norm_i_y, rho_Left, rho_Right, u_Left, u_Right, v_Left, v_Right,
		P_Left, P_Right, Fi_1, Fi_2, Fi_3, Fi_4);
      Calc_Flux(norm_j_x, norm_j_y, rho_Lower, rho_Upper, u_Lower, u_Upper, v_Lower, v_Upper,
		P_Lower, P_Upper, Fj_1, Fj_2, Fj_3, Fj_4);
      
      Calc_Res(Fi_1, Fj_1, S1, Ai, Aj, Volumn, Res1);
      Calc_Res(Fi_2, Fj_2, S2, Ai, Aj, Volumn, Res2);
      Calc_Res(Fi_3, Fj_3, S3, Ai, Aj, Volumn, Res3);
      Calc_Res(Fi_4, Fj_4, S4, Ai, Aj, Volumn, Res4);
    }
    //-------------------------------- * Runge-Kutta time marching END * ------------------------------------//
    
    // Calc L2 norms of the residuals
    L2norm1 = Calc_L2Norm(Res1);
    L2norm2 = Calc_L2Norm(Res2);
    L2norm3 = Calc_L2Norm(Res3);
    L2norm4 = Calc_L2Norm(Res4);

    L2norm1_rel.push_back(L2norm1/L2norm1_initial);
    L2norm2_rel.push_back(L2norm2/L2norm2_initial);
    L2norm3_rel.push_back(L2norm3/L2norm3_initial);
    L2norm4_rel.push_back(L2norm4/L2norm4_initial);

    if (flag_case ==1) {
      Calc_Error(rho, rho_MMS, Error1);
      Calc_Error(u, u_MMS, Error2);
      Calc_Error(v, v_MMS, Error3);
      Calc_Error(P, P_MMS, Error4);
      
      iterE_rho = Calc_L2Norm(Error1);
      iterE_u = Calc_L2Norm(Error2);
      iterE_v = Calc_L2Norm(Error3);
      iterE_P = Calc_L2Norm(Error4);

      /*
      iterError_rho.push_back(iterE_rho/iterE_rho_initial);
      iterError_u.push_back(iterE_u/iterE_u_initial);
      iterError_v.push_back(iterE_v/iterE_v_initial);
      iterError_P.push_back(iterE_P/iterE_P_initial);
      */

      iterError_rho.push_back(iterE_rho);
      iterError_u.push_back(iterE_u);
      iterError_v.push_back(iterE_v);
      iterError_P.push_back(iterE_P);
    }
   
    // save result to U_Rk to get ready for the next time step
    GoTomb_2d(0, U1_RK, U1);
    GoTomb_2d(0, U2_RK, U2);
    GoTomb_2d(0, U3_RK, U3);
    GoTomb_2d(0, U4_RK, U4);
    //---------------------------------------- Save data to file -------------------------------------------//
    
    if (((1.0*iteration/WriteInterval)-int(1.0*iteration/WriteInterval) ==0) ||
	 (iteration == 1 && debug == "true" ) )
      {
	Output_2D_Vector("rho_tomb", rho_tomb, iteration);
	Output_2D_Vector("u_tomb", u_tomb, iteration);
	Output_2D_Vector("v_tomb", v_tomb, iteration);
	Output_2D_Vector("P_tomb", P_tomb, iteration);

	Output_2D_Vector("rho", rho, iteration);
	Output_2D_Vector("u", u, iteration);
	Output_2D_Vector("v", v, iteration);
	Output_2D_Vector("P", P, iteration); 

	Output_2D_Vector("U1", U1, iteration);
	Output_2D_Vector("U2", U2, iteration);
	Output_2D_Vector("U3", U3, iteration);
	Output_2D_Vector("U4", U4, iteration);

	Output_2D_Vector("Fi_1", Fi_1, iteration);
	Output_2D_Vector("Fi_2", Fi_2, iteration);
	Output_2D_Vector("Fi_3", Fi_3, iteration);
	Output_2D_Vector("Fi_4", Fi_4, iteration);

	Output_2D_Vector("Fj_1", Fj_1, iteration);
	Output_2D_Vector("Fj_2", Fj_2, iteration);
	Output_2D_Vector("Fj_3", Fj_3, iteration);
	Output_2D_Vector("Fj_4", Fj_4, iteration);
	
	Output_2D_Vector("S1", S1, iteration);
	Output_2D_Vector("S2", S2, iteration);
	Output_2D_Vector("S3", S3, iteration);
	Output_2D_Vector("S4", S4, iteration);

	Output_2D_Vector("dt", dt, iteration);

	Output_2D_Vector("Res1", Res1, iteration);
	Output_2D_Vector("Res2", Res2, iteration);
	Output_2D_Vector("Res3", Res3, iteration);
	Output_2D_Vector("Res4", Res4, iteration);

	Output_2D_Vector("Error1", Error1, iteration);
	Output_2D_Vector("Error2", Error2, iteration);
	Output_2D_Vector("Error3", Error3, iteration);
	Output_2D_Vector("Error4", Error4, iteration);
	
	//===== Debug =====//
	if (true) {
	  Output_2D_Vector("rho_Left", rho_Left, iteration);
	  Output_2D_Vector("u_Left", u_Left, iteration);
	  Output_2D_Vector("v_Left", v_Left, iteration);
	  Output_2D_Vector("P_Left", P_Left, iteration);
	  Output_2D_Vector("rho_Right", rho_Right, iteration);
	  Output_2D_Vector("u_Right", u_Right, iteration);
	  Output_2D_Vector("v_Right", v_Right, iteration);
	  Output_2D_Vector("P_Right", P_Right, iteration);
	  Output_2D_Vector("rho_Lower", rho_Lower, iteration);
	  Output_2D_Vector("u_Lower", u_Lower, iteration);
	  Output_2D_Vector("v_Lower", v_Lower, iteration);
	  Output_2D_Vector("P_Lower", P_Lower, iteration);
	  Output_2D_Vector("rho_Upper", rho_Upper, iteration);
	  Output_2D_Vector("u_Upper", u_Upper, iteration);
	  Output_2D_Vector("v_Upper", v_Upper, iteration);
	  Output_2D_Vector("P_Upper", P_Upper, iteration);
	}	
      }
        
    //----------------------------------- Break if reaches the final residual -------------------------------//
    if ((L2norm1_rel.back() < Final_res) && (L2norm2_rel.back() < Final_res) &&
	(L2norm3_rel.back() < Final_res) && (L2norm4_rel.back() < Final_res))
      {
	Output_2D_Vector("rho_tomb", rho_tomb, 9999);
	Output_2D_Vector("u_tomb", u_tomb, 9999);
	Output_2D_Vector("v_tomb", v_tomb, 9999);
	Output_2D_Vector("P_tomb", P_tomb, 9999);

	Output_2D_Vector("rho", rho, 9999);
	Output_2D_Vector("u", u, 9999);
	Output_2D_Vector("v", v, 9999);
	Output_2D_Vector("P", P, 9999); 

	Output_2D_Vector("U1", U1, 9999);
	Output_2D_Vector("U2", U2, 9999);
	Output_2D_Vector("U3", U3, 9999);
	Output_2D_Vector("U4", U4, 9999);

	Output_2D_Vector("Fi_1", Fi_1, 9999);
	Output_2D_Vector("Fi_2", Fi_2, 9999);
	Output_2D_Vector("Fi_3", Fi_3, 9999);
	Output_2D_Vector("Fi_4", Fi_4, 9999);

	Output_2D_Vector("Fj_1", Fj_1, 9999);
	Output_2D_Vector("Fj_2", Fj_2, 9999);
	Output_2D_Vector("Fj_3", Fj_3, 9999);
	Output_2D_Vector("Fj_4", Fj_4, 9999);
	
	Output_2D_Vector("S1", S1, 9999);
	Output_2D_Vector("S2", S2, 9999);
	Output_2D_Vector("S3", S3, 9999);
	Output_2D_Vector("S4", S4, 9999);

	Output_2D_Vector("dt", dt, 9999);

	Output_2D_Vector("Res1", Res1, 9999);
	Output_2D_Vector("Res2", Res2, 9999);
	Output_2D_Vector("Res3", Res3, 9999);
	Output_2D_Vector("Res4", Res4, 9999);

	Output_2D_Vector("Error1", Error1, 9999);
	Output_2D_Vector("Error2", Error2, 9999);
	Output_2D_Vector("Error3", Error3, 9999);
	Output_2D_Vector("Error4", Error4, 9999);

	CoutTitle("Residuals", "~", "short");
	cout<<" * iteration:            "<<iteration<<endl;
	cout<<" * dt:                   "<<dt_global<<endl;
	cout<<" * Mass Equ Residual:    "<<L2norm1_rel.back()<<endl;
	cout<<" * x mmt Equ Residual:   "<<L2norm2_rel.back()<<endl;
	cout<<" * y mmt Equ Residual:   "<<L2norm3_rel.back()<<endl;
	cout<<" * Energy Equ Residual:  "<<L2norm4_rel.back()<<endl;
	break;
      }

    if ((1.0*iteration/PrintInterval)-int(1.0*iteration/PrintInterval) ==0) {
      CoutTitle("Input Variables", "~", "short");
      cout<<" * iteration:            "<<iteration<<endl;
      cout<<" * dt:                   "<<dt_global<<endl;
      cout<<" * Mass Equ Residual:    "<<L2norm1_rel.back()<<endl;
      cout<<" * x mmt Equ Residual:   "<<L2norm2_rel.back()<<endl;
      cout<<" * y mmt Equ Residual:   "<<L2norm3_rel.back()<<endl;
      cout<<" * Energy Equ Residual:  "<<L2norm4_rel.back()<<endl;
    }

    // calculate dt_global for the next time step
    Calc_Vs(rho, P, Vs);
    Calc_dt(CFL, Volumn, Ai, Aj, norm_i_x, norm_i_y, norm_j_x, norm_j_y, u, v, Vs, dt);
  }
  
  Output_1D_Vector("L2norm1_rel", L2norm1_rel);
  Output_1D_Vector("L2norm2_rel", L2norm2_rel);
  Output_1D_Vector("L2norm3_rel", L2norm3_rel);
  Output_1D_Vector("L2norm4_rel", L2norm4_rel);

  if (flag_case == 1) {
    Output_1D_Vector("iterError_rho", iterError_rho);
    Output_1D_Vector("iterError_u", iterError_u);
    Output_1D_Vector("iterError_v", iterError_v);
    Output_1D_Vector("iterError_P", iterError_P);
  }
  
  Save_data_to_folder();
}
  
//===============================================================================================================//
//                                          Other Numerical Functions                                            //
//===============================================================================================================//
void Calc_Area(int ni, int nj, vector<vector<double> > x, vector<vector<double> > y, vector<vector<double> > &Ai,
	       vector<vector<double> > &Aj)
{
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj-1; j++) {
      Ai[i][j] = sqrt( pow( (x[i][j] - x[i][j+1]), 2) + pow( (y[i][j] - y[i][j+1]),2) );
    }
  }
  for (int i=0; i<ni-1; i++) {
    for (int j=0; j<nj; j++) {
      Aj[i][j] = sqrt( pow( (x[i][j] - x[i+1][j]), 2) + pow( (y[i][j] - y[i+1][j]),2) );
    }
  }
      
}

void Calc_NormVec(vector<vector<double> > x, vector<vector<double> > y, vector<vector<double> > Ai,
		  vector<vector<double> > Aj, vector<vector<double> > &norm_i_x, vector<vector<double> > &norm_i_y,
		  vector<vector<double> > &norm_j_x, vector<vector<double> > &norm_j_y)
{
  int H_ni = Ai.size();
  int H_nj = Ai[0].size();
  int V_ni = Aj.size();
  int V_nj = Aj[0].size();

  for (int i=0; i<H_ni; i++) {
    for (int j=0; j<H_nj; j++) {
      norm_i_x[i][j] =  (y[i][j+1] - y[i][j]) / Ai[i][j];
      norm_i_y[i][j] = -(x[i][j+1] - x[i][j]) / Ai[i][j];
    }
  }

  for (int i=0; i<V_ni; i++) {
    for (int j=0; j<V_nj; j++) {
      norm_j_x[i][j] = -(y[i+1][j] - y[i][j]) / Aj[i][j];
      norm_j_y[i][j] = (x[i+1][j] - x[i][j]) / Aj[i][j];
    }
  }
}

void Calc_Volumn(int ni, int nj, vector<vector<double> > x, vector<vector<double> > y,
		 vector<vector<double> > &Volumn)
{
  double x1, y1, x2, y2;
  double cross_norm;
  for (int i=0; i<ni-1; i++) {
    for (int j=0; j<nj-1; j++) {
      x1 = x[i][j] - x[i+1][j+1];
      y1 = y[i][j] - y[i+1][j+1];
      x2 = x[i][j+1] - x[i+1][j];
      y2 = y[i][j+1] - y[i+1][j];
      cross_norm =  Vector_Cross_norm(x1, y1, x2, y2);
      Volumn[i][j] = 0.5*cross_norm;
    }
  }
}

void Resize_2dVec(int size_i, int size_j, vector<vector<double> > &Vec2d)
{
  Vec2d.resize(size_i);
  for (int i=0; i<size_i; i++)
    Vec2d[i].resize(size_j);
}

double Vector_Cross_norm(double x1, double y1, double x2, double y2)
{
  return abs(x1*y2-x2*y1);
}

//double Set_ghost_cell(vector<vector<double> >

void GoTomb_2d_extrapolate(int num_ghost, vector<vector<double> > Vec2d, vector<vector<double> > &Vec2d_tomb)
{
  int ni = Vec2d.size();
  int nj = Vec2d[0].size();
  for (int i=0; i<ni; i++) {
    copy(Vec2d[i].begin(), Vec2d[i].end(), Vec2d_tomb[i+num_ghost].begin()+num_ghost);
  }

  int lower_index = num_ghost - 1;
  int upper_index = num_ghost + nj;
  int left_index = num_ghost -1;
  int right_index = num_ghost + ni;

  for (int i=num_ghost; i<ni+num_ghost; i++) {
    for (int j=0; j<num_ghost; j++) {
      Vec2d_tomb[i][lower_index-j] = 2.0*Vec2d_tomb[i][lower_index-j+1] - Vec2d_tomb[i][lower_index-j+2];
      Vec2d_tomb[i][upper_index+j] = 2.0*Vec2d_tomb[i][upper_index+j-1] - Vec2d_tomb[i][upper_index+j-2];
    }
  }

  for (int j=num_ghost; j<nj+num_ghost; j++) {
    for (int i=0; i<num_ghost; i++) {
      Vec2d_tomb[left_index-i][j] = 2.0*Vec2d_tomb[left_index-i+1][j] - Vec2d_tomb[left_index-i+2][j];
      Vec2d_tomb[right_index+i][j] = 2.0*Vec2d_tomb[right_index+i-1][j] - Vec2d_tomb[right_index+i-2][j];
    }
  }
}

void GoTomb_2d(int num_ghost, vector<vector<double> > Vec2d, vector<vector<double> > &Vec2d_tomb)
{
  int ni = Vec2d.size();
  int nj = Vec2d[0].size();
  for (int i=0; i<ni; i++) {
    copy(Vec2d[i].begin(), Vec2d[i].end(), Vec2d_tomb[i+num_ghost].begin()+num_ghost);
  }
}

void OutTomb_2d(int num_ghost, vector<vector<double> > Vec2d_tomb, vector<vector<double> > &Vec2d)
{
  int ni = Vec2d.size();
  int nj = Vec2d[0].size();

  for (int i=0; i<ni; i++) {
    copy(Vec2d_tomb[i+num_ghost].begin() + num_ghost, Vec2d_tomb[i+num_ghost].end() - num_ghost, Vec2d[i].begin());
  }
}

void Ini_rho(vector<vector<double> > xc, vector<vector<double> > yc, vector<vector<double> > &rho)
{
  int ni = xc.size();
  int nj = xc[0].size();
  if (rho.size() != ni || rho[0].size() != nj)
    cout<<" Error: Ini_rho function, vectors size not match!!!"<<endl;

  double M;
  double rho0;
  double P0;
  double T0;
  double AOA;
  
  switch (flag_case)
    {
    case 1: //--------------------------------------- curviliniear --------------------------------------------//
      if (flag_supersonic == 0)        // subsonic
	{
	  rho0 = 1.0; 
	}
      else if (flag_supersonic == 1)   // supersonic
	{
	  rho0 = 1.0; 
	}
      else
	{
	  cout<<" *** Error: wrong input for flag_supersonic, choose 0 for subsonic or 1 for supersonic "<<endl;
	}
      for (int i=0; i<ni; i++) 
	for (int j=0; j<nj; j++) 
	  rho[i][j] = rho0;
      
      break;
    case 2: //-------------------------------------- 30 degree Inlet ------------------------------------------//
      M = 4.0;
      P0 = 12270.0;
      T0 = 217.0;
      rho0 = P0/(R*T0);
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  rho[i][j] = rho0;
      
      break;
    case 3: //---------------------------------------- NACA airfoil -------------------------------------------//
      if (flag_AOA == 1) {
	M = 0.84;  P0 = 65855.8;  T0 = 300.0;  AOA = 0;
      }
      if (flag_AOA == 2) {
	M = 0.82;  P0 = 67243.5;  T0 = 300.0;  AOA = 8;
	//M = 0.84;  P0 = 65855.8;  T0 = 300.0; AOA =8;
	
      }
      
      rho0 = P0/(R*T0);
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  rho[i][j] = rho0;
      
      break;
    } 
}

void Ini_u(vector<vector<double> > xc, vector<vector<double> > yc,
	   vector<vector<double> > &u)
{
  double uvel0;
  int ni = xc.size();
  int nj = xc[0].size();
  if (u.size() != ni || u[0].size() != nj)
    cout<<" Error: Ini_u function, vectors size not match!!!"<<endl;

  double M;
  double rho0;
  double P0;
  double T0;
  double AOA;
  
  switch (flag_case)
    {
    case 1: //--------------------------------------- curviliniear --------------------------------------------//
      
      if (flag_supersonic == 0)        // subsonic
	{
	  uvel0 = 70.0;
	}
      else if (flag_supersonic == 1)   // supersonic
	{
	  uvel0 = 800.0;
	}
      else
	{
	  cout<<" *** Error: wrong input for flag_supersonic, choose 0 for subsonic or 1 for supersonic "<<endl;
	}
      for (int i=0; i<ni; i++) 
	for (int j=0; j<nj; j++) 
	  u[i][j] = uvel0;
      
      break;
    case 2: //-------------------------------------- 30 degree Inlet ------------------------------------------//
      M = 4.0;
      P0 = 12270.0;
      T0 = 217.0;
      uvel0 = M*sqrt(gamma*R*T0);
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  u[i][j] = uvel0;

      break;
    case 3: //---------------------------------------- NACA airfoil -------------------------------------------//
      if (flag_AOA == 1) {
	M = 0.84;  P0 = 65855.8;  T0 = 300.0;  AOA = 0;
      }
      if (flag_AOA == 2) {
	M = 0.82;  P0 = 67243.5;  T0 = 300.0;  AOA = 8;
	//M = 0.84;  P0 = 65855.8;  T0 = 300.0; AOA =8;
	
      }
      
      uvel0 = M*sqrt(gamma*R*T0)*cos(AOA*Pi/180.0);
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  u[i][j] = uvel0;
      
      break;
    }

}

void Ini_v(vector<vector<double> > xc, vector<vector<double> > yc,
	   vector<vector<double> > &v)
{
  double vvel0;
  int ni = xc.size();
  int nj = xc[0].size();
  if (v.size() != ni || v[0].size() != nj)
    cout<<" Error: Ini_v function, vectors size not match!!!"<<endl;

  double M;
  double rho0;
  double P0;
  double T0;
  double AOA;
  
  switch (flag_case)
    {
    case 1: //--------------------------------------- curviliniear --------------------------------------------//
      
      if (flag_supersonic == 0)        // subsonic
	{
	  vvel0 = 90.0;
	}
      else if (flag_supersonic == 1)   // supersonic
	{
	  vvel0 = 800.0;
	}
      else
	{
	  cout<<" *** Error: wrong input for flag_supersonic, choose 0 for subsonic or 1 for supersonic "<<endl;
	}
      
      for (int i=0; i<ni; i++) 
	for (int j=0; j<nj; j++) 
	  v[i][j] = vvel0;
      
      break;
    case 2: //-------------------------------------- 30 degree Inlet ------------------------------------------//
      M = 4.0;
      P0 = 12270.0;
      T0 = 217.0;
      vvel0 = 0.0;
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  v[i][j] = vvel0;
      
      break;
    case 3: //---------------------------------------- NACA airfoil -------------------------------------------//
      if (flag_AOA == 1) {
	M = 0.84;  P0 = 65855.8;  T0 = 300.0;  AOA = 0;
      }
      if (flag_AOA == 2) {
	M = 0.82;  P0 = 67243.5;  T0 = 300.0;  AOA = 8;
	//M = 0.84;  P0 = 65855.8;  T0 = 300.0; AOA =8;
      }
      
      vvel0 = M*sqrt(gamma*R*T0)*sin(AOA*Pi/180.0);
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  v[i][j] = vvel0;
      
      break;
    }
}

void Ini_P(vector<vector<double> > xc, vector<vector<double> > yc, vector<vector<double> > &P)
{
  int ni = xc.size();
  int nj = xc[0].size();
  if (P.size() != ni || P[0].size() != nj)
    cout<<" Error: Ini_P function, vectors size not match!!!"<<endl;

  double M;
  double rho0;
  double P0;
  double T0;
  double AOA;
  
  switch (flag_case)
    {
    case 1: //--------------------------------------- curviliniear --------------------------------------------//
      if (flag_supersonic == 0)        // subsonic
	{
	  P0 = 1.0e5;
	}
      else if (flag_supersonic == 1)   // supersonic
	{
	  P0 = 1.0e5;
	}
      else
	{
	  cout<<" *** Error: wrong input for flag_supersonic, choose 0 for subsonic or 1 for supersonic "<<endl;
	}

      for (int i=0; i<ni; i++) 
	for (int j=0; j<nj; j++) 
	  P[i][j] = P0;
      
      break;
    case 2: //-------------------------------------- 30 degree Inlet ------------------------------------------//
      M = 4.0;
      P0 = 12270.0;
      T0 = 217.0;
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  P[i][j] = P0;
      
      break;
    case 3: //---------------------------------------- NACA airfoil -------------------------------------------//
      if (flag_AOA == 1) {
	M = 0.84;  P0 = 65855.8;  T0 = 300.0;  AOA = 0;
      }
      if (flag_AOA == 2) {
	M = 0.82;  P0 = 67243.5;  T0 = 300.0;  AOA = 8;
	//M = 0.84;  P0 = 65855.8;  T0 = 300.0; AOA =8;
	
      }
      
      for (int i=0; i<ni; i++)
	for (int j=0; j<nj; j++)
	  P[i][j] = P0;
      
      break;
    } 
}

void Calc_Temperature(vector<vector<double> > rho, vector<vector<double> > P, vector<vector<double> > &T)
{
  int ni = rho.size();
  int nj = rho[0].size();

  for (int i=0; i<ni; i++) 
    for (int j=0; j<nj; j++) 
      T[i][j] = P[i][j]/(R*rho[i][j]);
}

void solve_MMS(vector<vector<double> > xc, vector<vector<double> > yc, vector<vector<double> > &rho,
	       vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &P)
{
  double rho0, rhox, rhoy;
  double uvel0, uvelx, uvely;
  double vvel0, vvelx, vvely;
  double press0, pressx, pressy;

  double x, y;
  double L = 1.0;
  //double Pi = acos(-1);
  int ni = xc.size();
  int nj = xc[0].size();
  
  if (flag_supersonic == 0)        // subsonic
    {
      rho0 = 1.0;       rhox = 0.15;       rhoy = -0.1;
      uvel0 = 70.0;     uvelx = 5.0;       uvely = -7.0;
      vvel0 = 90.0;     vvelx = -15.0;     vvely = 8.5;
      press0 = 1.0e5;   pressx = 0.2e5;    pressy = 0.5e5;
    }
  else if (flag_supersonic == 1)   // supersonic
    {
      rho0 = 1.0;       rhox = 0.15;       rhoy = -0.1;
      uvel0 = 800.0;    uvelx = 50.0;      uvely = -30.0;
      vvel0 = 800.0;    vvelx = -75.0;     vvely = 40.0;
      press0 = 1.0e5;   pressx = 0.2e5;    pressy = 0.5e5;      
    }
  else
    {
      cout<<" *** Error: wrong input for flag_supersonic, choose 0 for subsonic or 1 for supersonic "<<endl;
    }
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      x = xc[i][j];
      y = yc[i][j];

      rho[i][j] = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);

      u[i][j] = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));

      v[i][j] = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));

      P[i][j] = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
    }
  }
}

void solve_MMS_source(vector<vector<double> > xc, vector<vector<double> > yc,
		      vector<vector<double> > &S1, vector<vector<double> > &S2, vector<vector<double> > &S3,
		      vector<vector<double> > &S4)
{
  double rho0, rhox, rhoy;
  double uvel0, uvelx, uvely;
  double vvel0, vvelx, vvely;
  double wvel0, wvelx, wvely;
  double press0, pressx, pressy;

  double x, y;
  double L = 1.0;
  //double Pi = acos(-1);
  int ni = xc.size();
  int nj = xc[0].size();

  if (S1.size() != ni || S1[0].size() != nj)
    cout<<" Error: solve_MMS_source function, vectors size not match!!!"<<endl;
  
  if (flag_supersonic == 0)        // subsonic
    {
      rho0 = 1.0;       rhox = 0.15;       rhoy = -0.1;
      uvel0 = 70.0;     uvelx = 5.0;       uvely = -7.0;
      vvel0 = 90.0;     vvelx = -15.0;     vvely = 8.5;
      press0 = 1.0e5;   pressx = 0.2e5;    pressy = 0.5e5;
    }
  else if (flag_supersonic == 1)   // supersonic
    {
      rho0 = 1.0;       rhox = 0.15;       rhoy = -0.1;
      uvel0 = 800.0;    uvelx = 50.0;      uvely = -30.0;
      vvel0 = 800.0;    vvelx = -75.0;     vvely = 40.0;
      press0 = 1.0e5;   pressx = 0.2e5;    pressy = 0.5e5;      
    }
  else
    {
      cout<<" *** Error: wrong input for flag_supersonic, choose 0 for subsonic or 1 for supersonic "<<endl;
    }
    
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      x = xc[i][j];
      y = yc[i][j];	
      
      S1[i][j] = (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))/(2.*L) +
	         (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))/(3.*L) +
	         (Pi*rhox*cos((Pi*x)/L)*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L -
	         (Pi*rhoy*sin((Pi*y)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
		 vvely*sin((2*Pi*y)/(3.*L))))/(2.*L);

      S2[i][j] = (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
		 (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L +
	         (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
		 (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/(3.*L) +
	         (Pi*rhox*cos((Pi*x)/L)*pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2))/L -
	         (2*Pi*pressx*sin((2*Pi*x)/L))/L - (Pi*rhoy*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
		 uvelx*sin((3*Pi*x)/(2.*L)))*sin((Pi*y)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
		 vvely*sin((2*Pi*y)/(3.*L))))/ (2.*L) - (3*Pi*uvely*(rho0 + rhoy*cos((Pi*y)/(2.*L)) +
		 rhox*sin((Pi*x)/L))*sin((3*Pi*y)/(5.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
		 vvely*sin((2*Pi*y)/(3.*L))))/(5.*L);

      S3[i][j] = (Pi*pressy*cos((Pi*y)/L))/L - (Pi*vvelx*sin((Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) +
		  rhox*sin((Pi*x)/L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/(2.*L) +
	         (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
		 (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(2.*L) +
	         (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
		 (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(3.*L) +
	         (Pi*rhox*cos((Pi*x)/L)*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*
		 (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/L -
	         (Pi*rhoy*sin((Pi*y)/(2.*L))*pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
		 vvely*sin((2*Pi*y)/(3.*L)),2))/(2.*L);
      
      S4[i][j] = (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*
	         ((-2*Pi*pressx*sin((2*Pi*x)/L))/L + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*
		 ((-2*Pi*pressx*sin((2*Pi*x)/L))/((-1 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.*L)) +
		 rhox*sin((Pi*x)/L))) +	((3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
		 uvelx*sin((3*Pi*x)/(2.*L))))/L - (Pi*vvelx*sin((Pi*x)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
		 vvely*sin((2*Pi*y)/(3.*L))))/L)/2. - (Pi*rhox*cos((Pi*x)/L)*(press0 + pressx*cos((2*Pi*x)/L) +
		 pressy*sin((Pi*y)/L)))/((-1 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.*L)) +
		 rhox*sin((Pi*x)/L),2))) + (Pi*rhox*cos((Pi*x)/L)*((pow(wvel0,2) + pow(uvel0 +
		 uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((Pi*x)/(2.*L))
		 + vvely*sin((2*Pi*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/
		 ((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))))/L) +
	         (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L) +
		 (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*((pow(wvel0,2) +
		 pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) +
		 pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. +
		 (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/
		 (2.*L)) + rhox*sin((Pi*x)/L))))))/(2.*L) + (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(press0 +
		 pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L) + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/
		 L))*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) +
		 pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. +
		 (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/
		 (2.*L)) + rhox*sin((Pi*x)/L))))))/(3.*L) + (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/
		 (3.*L)))*((Pi*pressy*cos((Pi*y)/L))/L - (Pi*rhoy*sin((Pi*y)/(2.*L))*((pow(wvel0,2) +
		 pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) +
		 pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. +
		 (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/
		 (2.*L)) + rhox*sin((Pi*x)/L)))))/(2.*L) + (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
		 rhox*sin((Pi*x)/L))*((Pi*pressy*cos((Pi*y)/L))/((-1 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.*L)) +
		 rhox*sin((Pi*x)/L))) + ((-6*Pi*uvely*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/
		 (2.*L)))*sin((3*Pi*y)/(5.*L)))/(5.*L) + (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*
		 (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(3.*L))/2. +
		 (Pi*rhoy*sin((Pi*y)/(2.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L)))/
		 (2.*(-1 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L),2))));
      
    }
  }
}

void Enforce_BC_MMS(int num_ghost, vector<vector<double> > prm_tomb_MMS, vector<vector<double> > &prm_tomb)
{
  int ni_tomb = prm_tomb.size();
  int nj_tomb = prm_tomb[0].size();

  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;

  for (int i=num_ghost; i<ni+num_ghost; i++)
    {
      for (int j=0; j<num_ghost; j++)
	{
	  prm_tomb[i][j] = prm_tomb_MMS[i][j];
	  prm_tomb[i][j+nj+num_ghost] = prm_tomb_MMS[i][j+nj+num_ghost];
	}
    }
  
  for (int i=0; i<num_ghost; i++)
    {
      for (int j=num_ghost; j<nj+num_ghost; j++)
	{
	  prm_tomb[i][j] = prm_tomb_MMS[i][j];
	  prm_tomb[i+ni+num_ghost][j] = prm_tomb_MMS[i+ni+num_ghost][j];
	}
    }
}

void Enforce_BC(int num_ghost, vector<vector<double> > ni_x, vector<vector<double> > ni_y,
		vector<vector<double> > nj_x, vector<vector<double> > nj_y,  vector<vector<double> > &rho_tomb,
		vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb,
		vector<vector<double> > &P_tomb, vector<vector<double> > &T_tomb)
{
  int ni = rho_tomb.size() - 2*num_ghost;         // number of cells on i direction
  int nj = rho_tomb[0].size() - 2*num_ghost;      // number of cells on j direction
  switch (flag_case)
    {
    case 1: //--------------------------------------- curviliniear --------------------------------------------//
      {
	cout<<"ERROR: Wrong Input of flag_case for Enforce_BC function !!!!"<<endl;      
	break;
      }
    case 2: //-------------------------------------- 30 degree Inlet ------------------------------------------//
      {
	int joint;
	joint = 10*pow(2, flag_mesh);
	
	int LowerWall_Begin[2] = {0, 0};
	int LowerWall_End[2] = {ni-1, 0};
	int UpperWall_Begin[2] = {joint, nj-1};
	int UpperWall_End[2] = {ni-1, nj-1};
	int Inlet_Begin[2] = {0, nj-1};
	int Inlet_End[2] = {joint-1, nj-1};
	int Outlet_Begin[2] = {ni-1, 0};
	int Outlet_End[2] = {ni-1, nj-1};
	int SymtrWall_Begin[2] = {0,0};
	int SymtrWall_End[2] = {0, nj-1};
	
	Symmetric_BC(SymtrWall_Begin, SymtrWall_End, num_ghost, rho_tomb, u_tomb, v_tomb, P_tomb);
	Inlet_BC(Inlet_Begin, Inlet_End, num_ghost, rho_tomb, u_tomb, v_tomb, P_tomb);
	Outlet_BC(Outlet_Begin, Outlet_End, num_ghost, rho_tomb, u_tomb, v_tomb, P_tomb);
	SlipWall_BC(LowerWall_Begin, LowerWall_End, num_ghost, ni_x, ni_y, nj_x, nj_y,
		    rho_tomb, u_tomb, v_tomb, P_tomb, T_tomb);
	SlipWall_BC(UpperWall_Begin, UpperWall_End, num_ghost, ni_x, ni_y, nj_x, nj_y,
		    rho_tomb, u_tomb, v_tomb, P_tomb, T_tomb);
	
	break;
      }
    case 3: //---------------------------------------- NACA airfoil -------------------------------------------//
      {
	int Airfoil_start = 4*pow(2, flag_mesh);
	int Airfoil_end = 20*pow(2, flag_mesh) -1;
	
	int FarField1_Begin[2] = {0, 0} ;
	int FarField1_End[2] = {0, nj-1};
	int FarField2_Begin[2] = {0, nj-1};
	int FarField2_End[2] = {ni-1, nj-1};
	int FarField3_Begin[2] = {ni-1, 0};
	int FarField3_End[2] = {ni-1, nj-1};
	int PeriodicA_Begin[2] = {0, 0};
	int PeriodicA_End[2] = {Airfoil_start-1, 0};
	int PeriodicB_Begin[2] = {ni-1, 0};
	int PeriodicB_End[2] = {Airfoil_end+1, 0};
	int Air_Begin[2] = {Airfoil_start, 0};
	int Air_End[2] = {Airfoil_end, 0};
	
	FarField_BC(FarField1_Begin, FarField1_End, num_ghost, rho_tomb, u_tomb, v_tomb, P_tomb);
	FarField_BC(FarField2_Begin, FarField2_End, num_ghost, rho_tomb, u_tomb, v_tomb, P_tomb);
	FarField_BC(FarField3_Begin, FarField3_End, num_ghost, rho_tomb, u_tomb, v_tomb, P_tomb);
	
	SlipWall_BC(Air_Begin, Air_End, num_ghost, ni_x, ni_y, nj_x, nj_y, rho_tomb, u_tomb, v_tomb, P_tomb,
		    T_tomb);
	
	Periodic_BC(PeriodicA_Begin, PeriodicA_End, PeriodicB_Begin, PeriodicB_End, num_ghost,
		    rho_tomb, u_tomb, v_tomb, P_tomb);
	
	break;
      }
    }
}

void Symmetric_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
		  vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb)
{
  int ni_tomb = rho_tomb.size();
  int nj_tomb = rho_tomb[0].size();
  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;
  int i_tomb_Begin = Begin[0] + num_ghost;
  int j_tomb_Begin = Begin[1] + num_ghost;
  int i_tomb_End = End[0] + num_ghost;
  int j_tomb_End = End[1] + num_ghost;
  int I, J;
  
  for (int j=0; j<=End[1]-Begin[1]; j++) {
    for (int i=0; i<num_ghost; i++) {
      I = i_tomb_Begin-1;  // i index of the ghost cell nearest to the boundary
      J = j_tomb_Begin;    // j index of the first cell along the boundary
      rho_tomb[I-i][J+j] = rho_tomb[I+i+1][J+j];
      u_tomb[I-i][J+j] = u_tomb[I+i+1][J+j];
      v_tomb[I-i][J+j] = -v_tomb[I+i+1][J+j];
      P_tomb[I-i][J+j] = P_tomb[I+i+1][J+j];
    }
  }
}

void Inlet_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
	      vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb)
{
  int ni_tomb = rho_tomb.size();
  int nj_tomb = rho_tomb[0].size();
  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;
  int i_tomb_Begin = Begin[0] + num_ghost;
  int j_tomb_Begin = Begin[1] + num_ghost;
  int i_tomb_End = End[0] + num_ghost;
  int j_tomb_End = End[1] + num_ghost;
  int I, J;
  
  double M;
  double uvel0;
  double vvel0;
  double rho0;
  double P0;
  double T0;
  M = 4.0;
  P0 = 12270.0;
  T0 = 217.0;
  rho0 = P0/(R*T0);
  uvel0 = M*sqrt(gamma*R*T0);
  vvel0 = 0.0;

  for (int i=0; i<=End[0]-Begin[0]; i++) {
    for (int j=0; j<num_ghost; j++) {
      I = i_tomb_Begin;    // i index of the first cell along the boundary
      J = j_tomb_Begin+1;  // j index of the ghost cell nearest to the boundary
      rho_tomb[I+i][J+j] = rho0;
      u_tomb[I+i][J+j] = uvel0;
      v_tomb[I+i][J+j] = vvel0;
      P_tomb[I+i][J+j] = P0;
    }
  }   
}

void Outlet_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
	       vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb)
{
  int ni_tomb = rho_tomb.size();
  int nj_tomb = rho_tomb[0].size();
  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;
  int i_tomb_Begin = Begin[0] + num_ghost;
  int j_tomb_Begin = Begin[1] + num_ghost;
  int i_tomb_End = End[0] + num_ghost;
  int j_tomb_End = End[1] + num_ghost;
  int I, J;

  for (int j=0; j<=End[1]-Begin[1]; j++) {
    for (int i=0; i<num_ghost; i++) {
      I = i_tomb_Begin+1;  // i index of the ghost cell nearest to the boundary
      J = j_tomb_Begin;    // j index of the first cell along the boundary
      u_tomb[I+i][J+j] = 2.0*u_tomb[I+i-1][J+j] - u_tomb[I+i-2][J+j];
      v_tomb[I+i][J+j] = 2.0*v_tomb[I+i-1][J+j] - v_tomb[I+i-2][J+j];
      P_tomb[I+i][J+j] = 2.0*P_tomb[I+i-1][J+j] - P_tomb[I+i-2][J+j];
      rho_tomb[I+i][J+j] = 2.0*rho_tomb[I+i-1][J+j] - rho_tomb[I+i-2][J+j];
      /*T_tomb[I+i][J+j] = 2.0*T_tomb[I+i-1][J+j] - T_tomb[I+i-2][J+j];
      rho_tomb[I+i][J+j] = P_tomb[I+i][J+j]/(R*T_tomb[I+i][J+j]);
      */
    }
  }
}

void SlipWall_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > ni_x,
		 vector<vector<double> > ni_y, vector<vector<double> > nj_x, vector<vector<double> > nj_y,
		 vector<vector<double> > &rho_tomb, vector<vector<double> > &u_tomb,
		 vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb, vector<vector<double> > &T_tomb)
{
  int ni_tomb = rho_tomb.size();
  int nj_tomb = rho_tomb[0].size();
  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;
  int i_tomb_Begin = Begin[0] + num_ghost;
  int j_tomb_Begin = Begin[1] + num_ghost;
  int i_tomb_End = End[0] + num_ghost;
  int j_tomb_End = End[1] + num_ghost;
  int I, J;
  int sign;
  int BC_index;
  double NX, NY, Uvel, Vvel;
  double Temperature;

  if (i_tomb_Begin == i_tomb_End) { // if BC along j direction;
    
    if (i_tomb_Begin == num_ghost) { // if ghost cells on the smaller i direction of the boundary
      sign = -1;
      BC_index = Begin[0];
    }
    else if (i_tomb_Begin == ni + num_ghost -1) {  // if ghost cells on the larger i direction of the boundary
      sign = 1;
      BC_index = Begin[0] + 1;
    }
    else
      cout<<" Error 1: wrong index for SlipWall_BC function"<<endl;
    
    for (int j=0; j<=End[1]-Begin[1]; j++) { // loop all cells on this boundary
      for (int i=0; i<num_ghost; i++) {
	I = i_tomb_Begin + sign*1;   // i index of the ghost cell nearest to the boundary
	J = j_tomb_Begin;            // j index of the first cell along the boundary
	NX = ni_x[BC_index][Begin[1]+j];
	NY = ni_y[BC_index][Begin[1]+j];
	Uvel = u_tomb[i_tomb_Begin-sign*i][J+j];
	Vvel = v_tomb[i_tomb_Begin-sign*i][J+j];
	Temperature = T_tomb[i_tomb_Begin-sign*i][J+j];
	
	u_tomb[I+sign*i][J+j] = -NX*(Uvel*NX + Vvel*NY) - NY*(-Uvel*NY + Vvel*NX);
	v_tomb[I+sign*i][J+j] = -NY*(Uvel*NX + Vvel*NY) + NX*(-Uvel*NY + Vvel*NX);
	
	P_tomb[I+sign*i][J+j] = 2.0*P_tomb[I+sign*i-sign*1][J+j] - P_tomb[I+sign*i-sign*2][J+j];

	T_tomb[I+sign*i][J+j] = Temperature;

	rho_tomb[I+sign*i][J+j] = P_tomb[I+sign*i][J+j]/(R*T_tomb[I+sign*i][J+j]);
      }
    }
  }

  else if (j_tomb_Begin == j_tomb_End) { // if BC along i direction

    if (j_tomb_Begin == num_ghost) { // if ghost cells on the smaller j direction of the boundary
      sign = -1;
      BC_index = Begin[1];
    }
    else if (j_tomb_Begin == nj + num_ghost -1) {  // if ghost cells on the larger j direction of the boundary
      sign = 1;
      BC_index = Begin[1] + 1;
    }
    else
      cout<<" Error 2: wrong index for SlipWall_BC function"<<endl;

    for (int i=0; i<=End[0]-Begin[0]; i++) {  // loop all cells on this boundary
      for (int j=0; j<num_ghost; j++) {
	I = i_tomb_Begin;          // i index of the first cell along the boundary
	J = j_tomb_Begin + sign*1; // j index of the ghost cell nearest to the boundary
	NX = nj_x[Begin[0]+i][BC_index];
	NY = nj_y[Begin[0]+i][BC_index];
	Uvel = u_tomb[I+i][j_tomb_Begin-sign*j];
	Vvel = v_tomb[I+i][j_tomb_Begin-sign*j];
	Temperature = T_tomb[I+i][j_tomb_Begin-sign*j];

	u_tomb[I+i][J+sign*j] = -NX*(Uvel*NX + Vvel*NY) - NY*(-Uvel*NY + Vvel*NX);
	v_tomb[I+i][J+sign*j] = -NY*(Uvel*NX + Vvel*NY) + NX*(-Uvel*NY + Vvel*NX);

	P_tomb[I+i][J+sign*j] = 2.0*P_tomb[I+i][J+sign*j-sign*1] - P_tomb[I+i][J+sign*j-sign*2];
	
	T_tomb[I+i][J+sign*j] = Temperature;

	rho_tomb[I+i][J+sign*j] = P_tomb[I+i][J+sign*j]/(R*T_tomb[I+i][J+sign*j]);
      }
    }
  }
  else
    cout<<" Error 3: wrong index for SlipWall_BC function"<<endl;
}

void FarField_BC(int Begin[], int End[], int num_ghost, vector<vector<double> > &rho_tomb,
		 vector<vector<double> > &u_tomb, vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb)
{
  int ni_tomb = rho_tomb.size();
  int nj_tomb = rho_tomb[0].size();
  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;
  int i_tomb_Begin = Begin[0] + num_ghost;
  int j_tomb_Begin = Begin[1] + num_ghost;
  int i_tomb_End = End[0] + num_ghost;
  int j_tomb_End = End[1] + num_ghost;
  int I, J;
  int sign;
  double M, rho0, uvel0, vvel0, P0, T0, AOA;
  
  if (flag_case != 3)
    cout<<" Error: Only flag_case =3 will call function FarField_BC!!! "<<endl;
  if (flag_AOA == 1) {
    M = 0.84;  P0 = 65855.8;  T0 = 300.0;  AOA = 0;
  }
  else if (flag_AOA == 2) {
    M = 0.82;  P0 = 67243.5;  T0 = 300.0;  AOA = 8;
    //M = 0.84;  P0 = 65855.8;  T0 = 300.0; AOA =8;
    
  }
  else
    cout<<" Error: flag_AOA can only be 1 or 2!!!"<<endl;

  rho0 = P0/(R*T0);
  uvel0 = M*sqrt(gamma*R*T0)*cos(AOA*Pi/180.0);
  vvel0 = M*sqrt(gamma*R*T0)*sin(AOA*Pi/180.0);

  //-----------------------
  if (i_tomb_Begin == i_tomb_End) { // if BC along j direction;
    
    if (i_tomb_Begin == num_ghost) { // if ghost cells on the smaller i direction of the boundary
      sign = -1;
      //      BC_index = Begin[0];
    }
    else if (i_tomb_Begin == ni + num_ghost -1) {  // if ghost cells on the larger i direction of the boundary
      sign = 1;
      //      BC_index = Begin[0] + 1;
    }
    else
      cout<<" Error 1: wrong index for SlipWall_BC function"<<endl;
    
    for (int j=0; j<=End[1]-Begin[1]; j++) { // loop all cells on this boundary
      for (int i=0; i<num_ghost; i++) {
	I = i_tomb_Begin + sign*1;   // i index of the ghost cell nearest to the boundary
	J = j_tomb_Begin;            // j index of the first cell along the boundary

	rho_tomb[I+sign*i][J+j] = rho0;
	u_tomb[I+sign*i][J+j] = uvel0;
	v_tomb[I+sign*i][J+j] = vvel0;
	P_tomb[I+sign*i][J+j] = P0;
      }
    }
  }
  
  //-----------------------
  else if (j_tomb_Begin == j_tomb_End) { // if BC along i direction
    
    if (j_tomb_Begin == num_ghost) { // if ghost cells on the smaller j direction of the boundary
      sign = -1;
      //      BC_index = Begin[1];
    }
    else if (j_tomb_Begin == nj + num_ghost -1) {  // if ghost cells on the larger j direction of the boundary
      sign = 1;
      //      BC_index = Begin[1] + 1;
    }
    else
      cout<<" Error 2: wrong index for SlipWall_BC function"<<endl;

    for (int i=0; i<=End[0]-Begin[0]; i++) {  // loop all cells on this boundary
      for (int j=0; j<num_ghost; j++) {
	I = i_tomb_Begin;          // i index of the first cell along the boundary
	J = j_tomb_Begin + sign*1; // j index of the ghost cell nearest to the boundary

	rho_tomb[I+i][J+sign*j] = rho0;
	u_tomb[I+i][J+sign*j] = uvel0;
	v_tomb[I+i][J+sign*j] = vvel0;
	P_tomb[I+i][J+sign*j] = P0;
      }
    }
  }
}

void Periodic_BC(int A_Begin[], int A_End[], int B_Begin[], int B_End[], int num_ghost,
		 vector<vector<double> > &rho_tomb, vector<vector<double> > &u_tomb,
		 vector<vector<double> > &v_tomb, vector<vector<double> > &P_tomb)
{
  int ni_tomb = rho_tomb.size();
  int nj_tomb = rho_tomb[0].size();
  int ni = ni_tomb - 2*num_ghost;
  int nj = nj_tomb - 2*num_ghost;
  
  int i_tomb_A_Begin = A_Begin[0] + num_ghost;
  int j_tomb_A_Begin = A_Begin[1] + num_ghost;
  int i_tomb_A_End = A_End[0] + num_ghost;
  int j_tomb_A_End = A_End[1] + num_ghost;
  
  int i_tomb_B_Begin = B_Begin[0] + num_ghost;
  int j_tomb_B_Begin = B_Begin[1] + num_ghost;
  int i_tomb_B_End = B_End[0] + num_ghost;
  int j_tomb_B_End = B_End[1] + num_ghost;
  
  int I_A, J_A, I_B, J_B;
  int sign;
  int B_direction;

  int Length;

  if (A_Begin[0] > A_End[0] || A_Begin[1] > A_End[1])
    cout<<" Error: Periodic_BC, A_Begin should have lower index number than A_End!!!"<<endl;
  
  //-----------------------
  if (i_tomb_A_Begin == i_tomb_A_End) { // if BC along j direction;
    
    if (i_tomb_A_Begin == num_ghost) { // if ghost cells on the smaller i direction of the boundary
      sign = -1;
    }
    else if (i_tomb_A_Begin == ni + num_ghost -1) {  // if ghost cells on the larger i direction of the boundary
      sign = 1;
    }
    else
      cout<<" Error 1: wrong index for Periodic_BC function"<<endl;

    // to be implemented
    cout<<" Error2: Periodic_BC function: BC along j direction option is not implemented!!!"<<endl;
  }
  
  //-----------------------
  else if (j_tomb_A_Begin == j_tomb_A_End) { // if BC along i direction
    Length = abs(i_tomb_A_Begin - i_tomb_A_End) + 1;

    if (Length != abs(i_tomb_B_Begin - i_tomb_B_End) + 1 )
      cout<<" Error 3: wrong coordinates of B boundary input for Periodic_BC function"<<endl;
    
    if (j_tomb_A_Begin == num_ghost) { // if ghost cells on the smaller j direction of the boundary
      sign = -1;
    }
    else if (j_tomb_A_Begin == nj + num_ghost -1) {  // if ghost cells on the larger j direction of the boundary
      sign = 1;
    }
    else
      cout<<" Error 4: wrong index for Periodic_BC function"<<endl;

    if (i_tomb_B_End > i_tomb_B_Begin)
      B_direction = 1;
    else if (i_tomb_B_End < i_tomb_B_Begin)
      B_direction = -1;
    else
      cout<<" Error 5: wrong Input of B boundary for Periodic_BC function"<<endl;
    
    for (int i=0; i<Length; i++) {  // loop all cells on this boundary
      for (int j=0; j<num_ghost; j++) {
	I_A = i_tomb_A_Begin;          // i index of the first cell along the boundary A
	I_B = i_tomb_B_Begin;          // i index of the first cell along the boundary B
	J_A = j_tomb_A_Begin + sign*1;   // j index of the ghost cell nearest to the boundary A
	J_B = j_tomb_B_Begin + sign*1;   // j index of the ghost cell nearest to the boundary B

	rho_tomb[I_A+i][J_A+sign*j] = rho_tomb[I_B+B_direction*i][j_tomb_B_Begin-sign*j];
	u_tomb[I_A+i][J_A+sign*j] = u_tomb[I_B+B_direction*i][j_tomb_B_Begin-sign*j];
	v_tomb[I_A+i][J_A+sign*j] = v_tomb[I_B+B_direction*i][j_tomb_B_Begin-sign*j];
	P_tomb[I_A+i][J_A+sign*j] = P_tomb[I_B+B_direction*i][j_tomb_B_Begin-sign*j];

	rho_tomb[I_B+B_direction*i][J_B+sign*j] = rho_tomb[I_A+i][j_tomb_A_Begin-sign*j];
	u_tomb[I_B+B_direction*i][J_B+sign*j] = u_tomb[I_A+i][j_tomb_A_Begin-sign*j];
	v_tomb[I_B+B_direction*i][J_B+sign*j] = v_tomb[I_A+i][j_tomb_A_Begin-sign*j];
	P_tomb[I_B+B_direction*i][J_B+sign*j] = P_tomb[I_A+i][j_tomb_A_Begin-sign*j];
      }
    }
  }
}

void Switch_Prm2Consrv(vector<vector<double> > V1, vector<vector<double> > V2,
		       vector<vector<double> > V3, vector<vector<double> > V4, vector<vector<double> > &U1,
		       vector<vector<double> > &U2, vector<vector<double> > &U3, vector<vector<double> > &U4)
{
  int ni = V1.size();
  int nj = V1[0].size();

  if (U1.size() != ni || U1[0].size() != nj)
    cout<<" Error: Switch_Prm2Consrv function, vectors size not match!!!"<<endl;

  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      U1[i][j] = V1[i][j];
      U2[i][j] = V1[i][j]*V2[i][j];
      U3[i][j] = V1[i][j]*V3[i][j];
      U4[i][j] = V4[i][j]/(gamma-1.0) + 0.5*( V1[i][j]*V2[i][j]*V2[i][j] + V1[i][j]*V3[i][j]*V3[i][j] );
    }
  }
}

void Switch_Consrv2Prm(vector<vector<double> > U1, vector<vector<double> > U2,
		       vector<vector<double> > U3, vector<vector<double> > U4, vector<vector<double> > &V1,
		       vector<vector<double> > &V2, vector<vector<double> > &V3, vector<vector<double> > &V4)
{
  int ni = U1.size();
  int nj = U1[0].size();

  if (V1.size() != ni || V1[0].size() != nj)
    cout<<" Error: Switch_Consrv2Prm function, vectors size not match!!!"<<endl;
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      V1[i][j] = U1[i][j];
      V2[i][j] = U2[i][j]/U1[i][j];
      V3[i][j] = U3[i][j]/U1[i][j];
      V4[i][j] = (gamma-1)*( U4[i][j] - 0.5*U2[i][j]*U2[i][j]/U1[i][j] - 0.5*U3[i][j]*U3[i][j]/U1[i][j] );
    }
  }
}

void Calc_upwind_Prm_i(int num_ghost, double upwind_eps, double upwind_kappa, vector<vector<double> > prm_tomb,
			vector<vector<double> > &prm_Left, vector<vector<double> > &prm_Right)
{
  double denom;
  double eps = 1e-6;
  double r1_pos, r2_pos, r2_neg, r3_neg;
  double phi1_pos, phi2_pos, phi2_neg, phi3_neg;

  int ni = prm_Left.size();
  int nj = prm_Left[0].size();
  
  int i_tomb, j_tomb;
  
  if (prm_tomb.size() != ni+2*num_ghost-1 || prm_tomb[0].size() != nj+2*num_ghost)
    cout<<" Error: Calc_upwind_Prim_i function, vectors size not match!!!"<<endl;

  for (int j=0; j<nj; j++) {
    j_tomb = j + num_ghost;
    for (int i=0; i<ni; i++) {
      i_tomb = i + num_ghost;    

      denom = (prm_tomb[i_tomb-1][j_tomb] - prm_tomb[i_tomb-2][j_tomb] );
      denom = SIGN( Max( abs(denom), eps ), denom );
      r1_pos = (prm_tomb[i_tomb][j_tomb] - prm_tomb[i_tomb-1][j_tomb])/denom;      // r[i-1/2]_pos

      denom = ( prm_tomb[i_tomb][j_tomb] - prm_tomb[i_tomb-1][j_tomb] );
      denom = SIGN( Max( abs(denom), eps ), denom );
      r2_pos = (prm_tomb[i_tomb+1][j_tomb] - prm_tomb[i_tomb][j_tomb])/denom;      // r[i+1/2]_pos
      r2_neg = (prm_tomb[i_tomb-1][j_tomb] - prm_tomb[i_tomb-2][j_tomb])/denom;    // r[i+1/2]_neg

      denom = ( prm_tomb[i_tomb+1][j_tomb] - prm_tomb[i_tomb][j_tomb] );
      denom = SIGN( Max( abs(denom), eps ), denom );
      r3_neg = (prm_tomb[i_tomb][j_tomb] - prm_tomb[i_tomb-1][j_tomb])/denom;      // r[i+3/2]_neg

      switch (flag_limiter)
	{
	case 0: // turn off Limiters
	  {
	    phi1_pos = 1.0;
	    phi2_pos = 1.0;
	    phi2_neg = 1.0;
	    phi3_neg = 1.0;
	  }
	case 1: // using Van Albada Limiters
	  {
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
	    break;
	  } 
	}
	  prm_Left[i][j] = prm_tomb[i_tomb-1][j_tomb] + 0.25*upwind_eps*((1-upwind_kappa) * phi1_pos *
	                   (prm_tomb[i_tomb-1][j_tomb]-prm_tomb[i_tomb-2][j_tomb]) +
		           (1+upwind_kappa) * phi2_neg * (prm_tomb[i_tomb][j_tomb]-prm_tomb[i_tomb-1][j_tomb]));
      
	  prm_Right[i][j] = prm_tomb[i_tomb][j_tomb] - 0.25*upwind_eps*((1+upwind_kappa) * phi2_pos *
			    (prm_tomb[i_tomb][j_tomb]-prm_tomb[i_tomb-1][j_tomb]) + (1-upwind_kappa) * phi3_neg *
			    (prm_tomb[i_tomb+1][j_tomb]-prm_tomb[i_tomb][j_tomb]));
      
    }
  }
}

void Calc_upwind_Prm_j(int num_ghost, double upwind_eps, double upwind_kappa, vector<vector<double> > prm_tomb,
			vector<vector<double> > &prm_Lower, vector<vector<double> > &prm_Upper)
{
  double denom;
  double eps = 1e-6;
  double r1_pos, r2_pos, r2_neg, r3_neg;
  double phi1_pos, phi2_pos, phi2_neg, phi3_neg;

  int ni = prm_Lower.size();
  int nj = prm_Lower[0].size();
  
  int i_tomb, j_tomb;
  
  if (prm_tomb.size() != ni+2*num_ghost || prm_tomb[0].size() != nj+2*num_ghost-1)
    cout<<" Error: Calc_upwind_Prim_j function, vectors size not match!!!"<<endl;
  
  for (int i=0; i<ni; i++) {
    i_tomb = i + num_ghost;    
    for (int j=0; j<nj; j++) {
      j_tomb = j + num_ghost;

      denom = ( prm_tomb[i_tomb][j_tomb-1] - prm_tomb[i_tomb][j_tomb-2] );
      denom = SIGN( Max( abs(denom), eps ), denom );
      r1_pos = (prm_tomb[i_tomb][j_tomb] - prm_tomb[i_tomb][j_tomb-1])/denom;     // r[i-1/2]_pos

      denom = ( prm_tomb[i_tomb][j_tomb] - prm_tomb[i_tomb][j_tomb-1] );
      denom = SIGN( Max( abs(denom), eps ), denom );
      r2_pos = (prm_tomb[i_tomb][j_tomb+1] - prm_tomb[i_tomb][j_tomb])/denom;     // r[i+1/2]_pos
      r2_neg = (prm_tomb[i_tomb][j_tomb-1] - prm_tomb[i_tomb][j_tomb-2])/denom;   // r[i+1/2]_neg

      denom = ( prm_tomb[i_tomb][j_tomb+1] - prm_tomb[i_tomb][j_tomb] );
      denom = SIGN( Max( abs(denom), eps ), denom );
      r3_neg = (prm_tomb[i_tomb][j_tomb] - prm_tomb[i_tomb][j_tomb-1])/denom;     // r[i+3/2]_neg

      switch (flag_limiter)
	{
	case 0: // turn off Limiters
	  {
	    phi1_pos = 1.0;
	    phi2_pos = 1.0;
	    phi2_neg = 1.0;
	    phi3_neg = 1.0;
	  }
	case 1: // using Van Albada Limiters
	  {  
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
	  }
	}
      
      prm_Lower[i][j] = prm_tomb[i_tomb][j_tomb-1] + 0.25*upwind_eps*((1-upwind_kappa) * phi1_pos *
		        (prm_tomb[i_tomb][j_tomb-1]-prm_tomb[i_tomb][j_tomb-2]) + (1+upwind_kappa) *
		        phi2_neg * (prm_tomb[i_tomb][j_tomb]-prm_tomb[i_tomb][j_tomb-1]));
      prm_Upper[i][j] = prm_tomb[i_tomb][j_tomb] - 0.25*upwind_eps*((1+upwind_kappa) * phi2_pos *
			(prm_tomb[i_tomb][j_tomb]-prm_tomb[i_tomb][j_tomb-1]) + (1-upwind_kappa) *
			phi3_neg * (prm_tomb[i_tomb][j_tomb+1]-prm_tomb[i_tomb][j_tomb]));
    }
  }
}

void MUSCL_prm(int num_ghost, double upwind_eps, double upwind_kappa, vector<vector<double> > rho_tomb,
	       vector<vector<double> > u_tomb, vector<vector<double> > v_tomb, vector<vector<double> > P_tomb,
	       vector<vector<double> > &rho_Left, vector<vector<double> > &rho_Right,
	       vector<vector<double> > &rho_Lower, vector<vector<double> > &rho_Upper,
	       vector<vector<double> > &u_Left, vector<vector<double> > &u_Right,
	       vector<vector<double> > &u_Lower, vector<vector<double> > &u_Upper,
	       vector<vector<double> > &v_Left, vector<vector<double> > &v_Right,
	       vector<vector<double> > &v_Lower, vector<vector<double> > &v_Upper,
	       vector<vector<double> > &P_Left, vector<vector<double> > &P_Right,
	       vector<vector<double> > &P_Lower, vector<vector<double> > &P_Upper)
{
  Calc_upwind_Prm_i(num_ghost, upwind_eps, upwind_kappa, rho_tomb, rho_Left, rho_Right);
  Calc_upwind_Prm_i(num_ghost, upwind_eps, upwind_kappa, u_tomb, u_Left, u_Right);
  Calc_upwind_Prm_i(num_ghost, upwind_eps, upwind_kappa, v_tomb, v_Left, v_Right);
  Calc_upwind_Prm_i(num_ghost, upwind_eps, upwind_kappa, P_tomb, P_Left, P_Right);

  Calc_upwind_Prm_j(num_ghost, upwind_eps, upwind_kappa, rho_tomb, rho_Lower, rho_Upper);
  Calc_upwind_Prm_j(num_ghost, upwind_eps, upwind_kappa, u_tomb, u_Lower, u_Upper);
  Calc_upwind_Prm_j(num_ghost, upwind_eps, upwind_kappa, v_tomb, v_Lower, v_Upper);
  Calc_upwind_Prm_j(num_ghost, upwind_eps, upwind_kappa, P_tomb, P_Lower, P_Upper);
}

void Calc_Flux_VL(vector<vector<double> > nx, vector<vector<double> > ny,
		  vector<vector<double> > rho_L, vector<vector<double> > rho_R, vector<vector<double> > u_L,
		  vector<vector<double> > u_R, vector<vector<double> > v_L, vector<vector<double> > v_R,
		  vector<vector<double> > P_L, vector<vector<double> > P_R, vector<vector<double> > &F1,
		  vector<vector<double> > &F2, vector<vector<double> > &F3, vector<vector<double> > &F4)
{
  double M_L, M_R, M_pos, M_neg;
  double beta_L, beta_R;
  double alpha_pos, alpha_neg;
  double c_pos, c_neg;
  double Vs_L, Vs_R;
  double U_hec_L, U_hec_R; 
  double D_pos, D_neg;
  double Pbar_pos, Pbar_neg;
  double ht_L, ht_R;
  double Fc1, Fc2, Fc3, Fc4;
  double Fp1, Fp2, Fp3, Fp4;

  int ni = rho_L.size();
  int nj = rho_L[0].size();

  if (nx.size() != ni || nx[0].size() != nj)
    cout<<" Error: Calc_Flux_VL function, vectors size not match!!!"<<endl;
  
  for (int i = 0; i<ni; i++) {
    for (int j = 0; j<nj; j++) {
      Vs_L = sqrt(gamma*P_L[i][j]/rho_L[i][j]);
      Vs_R = sqrt(gamma*P_R[i][j]/rho_R[i][j]);

      U_hec_L = u_L[i][j]*nx[i][j] + v_L[i][j]*ny[i][j];
      U_hec_R = u_R[i][j]*nx[i][j] + v_R[i][j]*ny[i][j];

      M_L = U_hec_L/Vs_L;
      M_R = U_hec_R/Vs_R;

      M_pos = 0.25*(M_L + 1)*(M_L + 1);
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
      
      ht_L = (gamma/(gamma-1))*P_L[i][j]/rho_L[i][j] + 0.5*( u_L[i][j]*u_L[i][j] + v_L[i][j]*v_L[i][j] );
      ht_R = (gamma/(gamma-1))*P_R[i][j]/rho_R[i][j] + 0.5*(u_R[i][j]*u_R[i][j] + v_R[i][j]*v_R[i][j] );

      Fc1 = rho_L[i][j]*Vs_L*c_pos + rho_R[i][j]*Vs_R*c_neg;
      Fc2 = rho_L[i][j]*Vs_L*c_pos*u_L[i][j] + rho_R[i][j]*Vs_R*c_neg*u_R[i][j];
      Fc3 = rho_L[i][j]*Vs_L*c_pos*v_L[i][j] + rho_R[i][j]*Vs_R*c_neg*v_R[i][j];
      Fc4 = rho_L[i][j]*Vs_L*c_pos*ht_L + rho_R[i][j]*Vs_R*c_neg*ht_R;

      Fp1 = 0;
      Fp2 = D_pos*nx[i][j]*P_L[i][j] + D_neg*nx[i][j]*P_R[i][j];
      Fp3 = D_pos*ny[i][j]*P_L[i][j] + D_neg*ny[i][j]*P_R[i][j];
      Fp4 = 0;

      F1[i][j] = Fc1 + Fp1;
      F2[i][j] = Fc2 + Fp2;
      F3[i][j] = Fc3 + Fp3;
      F4[i][j] = Fc4 + Fp4;
    }
  } 
}

void Calc_Flux_Roe(vector<vector<double> > nx, vector<vector<double> > ny,
		   vector<vector<double> > rho_L, vector<vector<double> > rho_R, vector<vector<double> > u_L,
		   vector<vector<double> > u_R, vector<vector<double> > v_L, vector<vector<double> > v_R,
		   vector<vector<double> > P_L, vector<vector<double> > P_R, vector<vector<double> > &F1,
		   vector<vector<double> > &F2, vector<vector<double> > &F3, vector<vector<double> > &F4)
{
  int ne = 4;           // Number of equations we have in our Euler system

  double eps = 0.1;     // the eps for the modified Roe averaged eigenvalues
  double rho_Roe;
  double u_Roe;
  double v_Roe;
  double U_hat_Roe_L, U_hat_Roe_R;
  double U_hat_Roe;
  double ht_L, ht_R;
  double ht_Roe;
  double Vs_Roe;        // Roe averaged sound speed
  double R_Roe;         // sqrt(rho_R/rho_L)
  double drho, du, dv, dP;

  double abs_lambda1_Roe, abs_lambda2_Roe, abs_lambda3_Roe, abs_lambda4_Roe;
  double r1_Roe[ne], r2_Roe[ne], r3_Roe[ne], r4_Roe[ne];
  double dw1, dw2, dw3, dw4;

  double F1_L, F1_R, F2_L, F2_R, F3_L, F3_R, F4_L, F4_R;
  double sum[ne];
  
  int ni = rho_L.size();
  int nj = rho_L[0].size();

  if (nx.size() != ni || nx[0].size() != nj)
    cout<<" Error: Calc_Flux_VL function, vectors size not match!!!"<<endl;
  
  for (int i = 0; i<ni; i++) {
    for (int j = 0; j<nj; j++) {
      R_Roe = sqrt(rho_R[i][j]/rho_L[i][j]);
      rho_Roe = R_Roe*rho_L[i][j];
      u_Roe = (R_Roe*u_R[i][j] + u_L[i][j]) / (R_Roe + 1);
      v_Roe = (R_Roe*v_R[i][j] + v_L[i][j]) / (R_Roe + 1);
      U_hat_Roe = u_Roe*nx[i][j] + v_Roe*ny[i][j];
      
      ht_L = (gamma/(gamma-1))*P_L[i][j]/rho_L[i][j] + 0.5*( u_L[i][j]*u_L[i][j] + v_L[i][j]*v_L[i][j] );
      ht_R = (gamma/(gamma-1))*P_R[i][j]/rho_R[i][j] + 0.5*( u_R[i][j]*u_R[i][j] + v_R[i][j]*v_R[i][j] );
      ht_Roe = (R_Roe*ht_R + ht_L) / (R_Roe + 1);

      Vs_Roe = sqrt( (gamma-1)*( ht_Roe - 0.5*(u_Roe*u_Roe + v_Roe*v_Roe) ) );

      abs_lambda1_Roe = abs(U_hat_Roe);
      abs_lambda2_Roe = abs(U_hat_Roe);
      abs_lambda3_Roe = abs(U_hat_Roe + Vs_Roe);
      abs_lambda4_Roe = abs(U_hat_Roe - Vs_Roe);
      
      // calc modified lambda
      if (abs_lambda1_Roe < 2*eps*Vs_Roe) {
	abs_lambda1_Roe = abs_lambda1_Roe*abs_lambda1_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;
	abs_lambda2_Roe = abs_lambda2_Roe*abs_lambda2_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;
      }
      
      if (abs_lambda3_Roe < 2*eps*Vs_Roe)
	abs_lambda3_Roe = abs_lambda3_Roe*abs_lambda3_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;
      if (abs_lambda4_Roe < 2*eps*Vs_Roe)
	abs_lambda4_Roe = abs_lambda4_Roe*abs_lambda4_Roe/(4*eps*Vs_Roe) + eps*Vs_Roe;
      
      r1_Roe[0] = 1.0;
      r1_Roe[1] = u_Roe;
      r1_Roe[2] = v_Roe;
      r1_Roe[3] = 0.5*( u_Roe*u_Roe + v_Roe*v_Roe  );

      r2_Roe[0] = 0;
      r2_Roe[1] = rho_Roe * ny[i][j];
      r2_Roe[2] = -rho_Roe * nx[i][j];
      r2_Roe[3] = rho_Roe*( u_Roe*ny[i][j] - v_Roe*nx[i][j] );
	
      r3_Roe[0] = (0.5*rho_Roe/Vs_Roe);
      r3_Roe[1] = (0.5*rho_Roe/Vs_Roe) * ( u_Roe + Vs_Roe*nx[i][j] );
      r3_Roe[2] = (0.5*rho_Roe/Vs_Roe) * ( v_Roe + Vs_Roe*ny[i][j] );
      r3_Roe[3] = (0.5*rho_Roe/Vs_Roe) * ( ht_Roe + Vs_Roe*U_hat_Roe );

      r4_Roe[0] = (-0.5*rho_Roe/Vs_Roe);
      r4_Roe[1] = (-0.5*rho_Roe/Vs_Roe) * ( u_Roe - Vs_Roe*nx[i][j] );
      r4_Roe[2] = (-0.5*rho_Roe/Vs_Roe) * ( v_Roe - Vs_Roe*ny[i][j] );
      r4_Roe[3] = (-0.5*rho_Roe/Vs_Roe) * ( ht_Roe - Vs_Roe*U_hat_Roe );

      drho = rho_R[i][j] - rho_L[i][j];
      du = u_R[i][j] - u_L[i][j];
      dv = v_R[i][j] - v_L[i][j];
      dP = P_R[i][j] - P_L[i][j];

      dw1 = drho - dP/(Vs_Roe*Vs_Roe);
      dw2 = du*ny[i][j] - dv*nx[i][j];
      dw3 = du*nx[i][j] + dv*ny[i][j] + dP/(rho_Roe*Vs_Roe);
      dw4 = du*nx[i][j] + dv*ny[i][j] - dP/(rho_Roe*Vs_Roe);

      U_hat_Roe_L = ( u_L[i][j]*nx[i][j] + v_L[i][j]*ny[i][j] );
      U_hat_Roe_R = ( u_R[i][j]*nx[i][j] + v_R[i][j]*ny[i][j] );
      
      F1_L = rho_L[i][j]*U_hat_Roe_L;
      F1_R = rho_R[i][j]*U_hat_Roe_R;

      F2_L = rho_L[i][j]*u_L[i][j]*U_hat_Roe_L + P_L[i][j]*nx[i][j];
      F2_R = rho_R[i][j]*u_R[i][j]*U_hat_Roe_R + P_R[i][j]*nx[i][j];

      F3_L = rho_L[i][j]*v_L[i][j]*U_hat_Roe_L + P_L[i][j]*ny[i][j];
      F3_R = rho_R[i][j]*v_R[i][j]*U_hat_Roe_R + P_R[i][j]*ny[i][j];

      F4_L = rho_L[i][j]*ht_L*U_hat_Roe_L;
      F4_R = rho_R[i][j]*ht_R*U_hat_Roe_R;

      for (int k=0; k<ne; k++) {
	sum[k] = 0.5*( abs(abs_lambda1_Roe)*dw1*r1_Roe[k] + abs(abs_lambda2_Roe)*dw2*r2_Roe[k]
		    + abs(abs_lambda3_Roe)*dw3*r3_Roe[k] + abs(abs_lambda4_Roe)*dw4*r4_Roe[k] );
      }
      
      F1[i][j] = 0.5*(F1_L + F1_R) - sum[0];
      F2[i][j] = 0.5*(F2_L + F2_R) - sum[1];
      F3[i][j] = 0.5*(F3_L + F3_R) - sum[2];
      F4[i][j] = 0.5*(F4_L + F4_R) - sum[3];
    }
  }	
}

void Calc_Flux(vector<vector<double> > nx, vector<vector<double> > ny,
	       vector<vector<double> > rho_L, vector<vector<double> > rho_R, vector<vector<double> > u_L,
	       vector<vector<double> > u_R, vector<vector<double> > v_L, vector<vector<double> > v_R,
	       vector<vector<double> > P_L, vector<vector<double> > P_R, vector<vector<double> > &F1,
	       vector<vector<double> > &F2, vector<vector<double> > &F3, vector<vector<double> > &F4)
{
  switch (flag_upwind)
    {
    case 1:                 // Van Leer
      Calc_Flux_VL(nx, ny, rho_L, rho_R, u_L, u_R, v_L, v_R, P_L, P_R, F1, F2, F3, F4);
      break;
      
    case 2:                 // Roe
      Calc_Flux_Roe(nx, ny, rho_L, rho_R, u_L, u_R, v_L, v_R, P_L, P_R, F1, F2, F3, F4);
      break;
    }
}

void Calc_Vs(vector<vector<double> > rho, vector<vector<double> > P, vector<vector<double> > &Vs)
{
  int ni = Vs.size();
  int nj = Vs[0].size();

  if (rho.size() != ni || rho[0].size() != nj)
    cout<<" Error: Calc_Vs function, vectors size not match!!!"<<endl;

  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      Vs[i][j] = sqrt(gamma*P[i][j]/rho[i][j]);
    }
  }
}

void Calc_dt(double CFL, vector<vector<double> > Volumn, vector<vector<double> > Ai,
	     vector<vector<double> > Aj, vector<vector<double> > ni_x, vector<vector<double> > ni_y,
	     vector<vector<double> > nj_x, vector<vector<double> > nj_y, vector<vector<double> > u,
	     vector<vector<double> > v, vector<vector<double> > Vs, vector<vector<double> > &dt)
{
  int ni = dt.size();
  int nj = dt[0].size();

  if (Volumn.size() != ni || Volumn[0].size() != nj)
    cout<<" Error: Calc_dt function, vectors size not match!!!"<<endl;
  
  double ni_hat_x, ni_hat_y, nj_hat_x, nj_hat_y;
  double lambda_i, lambda_j;
  double A_hat_i, A_hat_j;
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      ni_hat_x = 0.5*(ni_x[i][j] + ni_x[i+1][j]);
      ni_hat_y = 0.5*(ni_y[i][j] + ni_y[i+1][j]);
      nj_hat_x = 0.5*(nj_x[i][j] + nj_x[i][j+1]);
      nj_hat_y = 0.5*(nj_y[i][j] + nj_y[i][j+1]);

      A_hat_i = 0.5*(Ai[i][j] + Ai[i+1][j]);
      A_hat_j = 0.5*(Aj[i][j] + Aj[i][j+1]);
      
      lambda_i = abs( u[i][j]*ni_hat_x + v[i][j]*ni_hat_y ) + Vs[i][j];
      lambda_j = abs( u[i][j]*nj_hat_x + v[i][j]*nj_hat_y ) + Vs[i][j];

      dt[i][j] = CFL*Volumn[i][j]/( lambda_i*A_hat_i + lambda_j*A_hat_j );
    }
  }
}
 
void Calc_Res(vector<vector<double> > Fi, vector<vector<double> > Fj, vector<vector<double> > S,
	      vector<vector<double> > Ai, vector<vector<double> > Aj, vector<vector<double> > volumn,
	      vector<vector<double> > &Res)
{
  int ni = Res.size();
  int nj = Res[0].size();

  if (Fi.size() != ni+1 || Fi[0].size() != nj)
    cout<<" Error: Calc_Res function, vectors size not match!!!"<<endl;
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      Res[i][j] = Fi[i+1][j]*Ai[i+1][j] - Fi[i][j]*Ai[i][j] + Fj[i][j+1]*Aj[i][j+1] - Fj[i][j]*Aj[i][j]
	- S[i][j]*volumn[i][j];
    }
  }
}

void Calc_Error(vector<vector<double> > Simulation, vector<vector<double> > Exact, vector<vector<double> > &Error)
{
  int ni = Simulation.size();
  int nj = Simulation[0].size();

  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      Error[i][j] = Simulation[i][j] - Exact[i][j];
    }
  }
}
 
void RK_marching(int stage, vector<vector<double> > Consrv_old, double dt_global, vector<vector<double> > Volumn,
		 vector<vector<double> > Res, vector<vector<double> > &Consrv)
{
  int ni = Consrv.size();
  int nj = Consrv[0].size();

  if (Res.size() != ni || Res[0].size() != nj)
    cout<<" Error: Calc_dt function, vectors size not match!!!"<<endl;
    
  double a[4];

  if (RK_order == 1) {
    a[0] = 1.0; a[1]=a[2]=a[3]=0;
  }
  if (RK_order == 2) {
    a[0] = 0.5; a[1] = 1.0; a[2]=a[3]=0;
  }
  if (RK_order == 3) {
    cout<<"3rd Order Runge-Kutta scheme not available";
  }
  if (RK_order ==4) {
    a[0] = 0.25; a[1] = 1/3; a[2] = 0.5; a[3] =1.0;
  }
  if (RK_order >4) {
    cout<<"Runge-Kutta scheme not available for order higher than 4";
  }

  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      Consrv[i][j] = Consrv_old[i][j] - a[stage]*dt_global*Res[i][j]/Volumn[i][j];
    }
  }  
}

double Calc_L2Norm(vector<vector<double> > X)
{
  int ni = X.size();
  int nj = X[0].size();
  double sum = 0;

  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      sum = sum + X[i][j]*X[i][j];
    }
  }
  return sqrt(sum/(ni*nj));
}

double Max(double a, double b)
{
  return ( 0.5*(a+b + abs(a-b)) );
}

double Min(double a, double b)
{
  return ( 0.5*(a+b - abs(a-b)) );
}

double Min_Element(vector<vector<double> > Vec)
{
  int ni = Vec.size();
  int nj = Vec[0].size();
  double Min_Element = Vec[0][0];
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      Min_Element = Min(Min_Element, Vec[i][j]);
    }
  }
  return Min_Element;
}

double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}

//===============================================================================================================//
//                                                 I/O Functions                                                 //
//===============================================================================================================//
double load_inputNum(string FileName, string Var)
{
  const char *FILENAME = FileName.c_str();

  ifstream in(FILENAME);    
  stringstream buffer;
  buffer << in.rdbuf();
  string text = buffer.str();
  
  // define cursors
  size_t pos1 = 0;
  size_t pos2;
  
  //create the array to store the strings.
  string InputStr;
  string InputName = Var + ":";
  
  pos1 = text.find(InputName, pos1); //search for the test in " ". pos1 will be where the test was found.
  pos1 = text.find(":", pos1);
  pos2 = text.find(";", pos1); //search for the semicolon ";". pos2 will be where the test was found.
  InputStr = text.substr(pos1+1, (pos2-pos1-1));
  
  FILENAME = new char; delete FILENAME; FILENAME = NULL;
  
  return atof(InputStr.c_str()); 
}

string read_mesh_name()
{
  if (flag_case == 1) { // curve linear
    if (flag_mesh == 1)
      return "./Curviliniear_Grids/curv2d9.txt";
    if (flag_mesh == 2)
      return "./Curviliniear_Grids/curv2d17.txt";
    if (flag_mesh == 3)
      return "./Curviliniear_Grids/curv2d33.txt";
    if (flag_mesh == 4)
      return "./Curviliniear_Grids/curv2d65.txt";
    if (flag_mesh == 5)
      return "./Curviliniear_Grids/curv2d129.txt";
    if (flag_mesh == 6)
      return "./Curviliniear_Grids/curv2d257.txt";
  }
  if (flag_case == 2) { // 30 degree inlet
    if (flag_mesh == 1)
      return "./Inlet_Grids/Inlet.53x17.txt";
    if (flag_mesh == 2)
      return "./Inlet_Grids/Inlet.105x33.txt";
    if (flag_mesh == 3)
      return "./Inlet_Grids/Inlet.209x65.txt";
    if (flag_mesh == 4)
      return "./Inlet_Grids/Inlet.417x129.txt";
  }
  if (flag_case == 3) { // NACA airfoil
    if (flag_mesh == 1)
      return "./NACA64A006_Grids/NACA64A006.49x14.txt";
    if (flag_mesh == 2)
      return "./NACA64A006_Grids/NACA64A006.97x27.txt";
    if (flag_mesh == 3)
      return "./NACA64A006_Grids/NACA64A006.193x53.txt";
    if (flag_mesh == 4)
      return "./NACA64A006_Grids/NACA64A006.385x105.txt";
  }
  return "failure";
}
  
void InputMesh(int &ni, int &nj, int &nk, vector<vector<double> > &x, vector<vector<double> > &xc,
	       vector<vector<double> > &y, vector<vector<double> > &yc, vector<vector<double> > &z,
	       vector<vector<double> > &zc)
{
  string line;
  double value;
  vector<double> data;
  string filename = read_mesh_name();
  const char *FILENAME = filename.c_str();
  ifstream fileIN;
  //  int nblocks;
  
  fileIN.open(FILENAME);

  // Error Check
  if (fileIN.fail()) {
    cerr << "* File you are trying to access cannot be found or opened!\n";
    exit(1);
  }
  
  // Read data from file
  while (fileIN.good()) {
    while (getline(fileIN, line)) {
      istringstream streamA(line);
      while (streamA >> value) {
	data.push_back(value);
      }
    }
  }
  /*  for (int f=0; f<data.size(); f++)
      cout<<"data["<<f<<"] = "<<data[f]<<endl; */
  ni = data[1]-1; nj = data[2]-1; nk = data[3]-1;
  x.resize(ni+1); y.resize(ni+1); z.resize(ni+1);
  xc.resize(ni); yc.resize(ni); zc.resize(ni);
  for (int i=0; i<ni+1; i++) {
    x[i].resize(nj+1);
    y[i].resize(nj+1);
    z[i].resize(nj+1);
  }
  for (int i=0; i<ni; i++) {
    xc[i].resize(nj);
    yc[i].resize(nj);
    zc[i].resize(nj);
  }
  
  int x_index = 4;
  int y_index = x_index + (ni+1)*(nj+1)*(nk+1);
  int z_index = y_index + (ni+1)*(nj+1)*(nk+1);
  
  for (int j=0; j<nj+1; j++) {
    for (int i=0; i<ni+1; i++) {
      x[i][j] = data[x_index];
      //cout << x[i][j] << " i " << i << " j " << j << endl;
      y[i][j] = data[y_index];
      z[i][j] = data[z_index];
      x_index++;
      y_index++;
      z_index++;
    }
  }

  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      xc[i][j] = 0.25*(x[i][j] + x[i+1][j] + x[i][j+1] + x[i+1][j+1]);
      yc[i][j] = 0.25*(y[i][j] + y[i+1][j] + y[i][j+1] + y[i+1][j+1]);
      zc[i][j] = 0.25*(z[i][j] + z[i+1][j] + z[i][j+1] + z[i+1][j+1]);
    }
  }

  FILENAME = new char; delete FILENAME; FILENAME = NULL;
}

void build_case_folder()
{
  ostringstream StrConvert1;
  ostringstream StrConvert2;
  ostringstream StrConvert3;
  ostringstream StrConvert4;
  ostringstream StrConvert5;
  StrConvert1 << flag_case;
  string Case_index = StrConvert1.str();
  StrConvert2 << flag_mesh;
  string Mesh_index = StrConvert2.str();
  StrConvert3 << flag_upwind;
  string Flux_index = StrConvert3.str();
  StrConvert4 << upwind_eps+1;
  string Order_index = StrConvert4.str();
  StrConvert5 << flag_AOA;
  string AOA_index = StrConvert5.str();

  if (flag_case != 3)
    Output_Folder= "Case_" + Case_index + "_Mesh_" + Mesh_index + "_upwind_" + Flux_index + "_order_"+ Order_index;
  if (flag_case == 3)
    Output_Folder= "Case_" + Case_index + "_Mesh_" + Mesh_index + "_upwind_" + Flux_index + "_AOA_"+ AOA_index;
  string clear_folder = "exec rm -f ./"+Output_Folder+"/*";
  string make_folder = "exec mkdir -p " + Output_Folder;
  const char *execline1 = clear_folder.c_str();
  const char *execline2 = make_folder.c_str();
  system(execline1);
  system(execline2);

  execline1 = new char; delete execline1; execline1 = NULL;
  execline2 = new char; delete execline2; execline2 = NULL; 
}

void Output_2D_Vector(string FileName, vector<vector<double> > Vec2d, int iteration)
{
  ostringstream StrConvert;
  StrConvert << iteration;
  string num_iter = StrConvert.str();
  string dash = "-";
  string Address = "./" + Output_Folder + "/";
  string Suffix = ".txt";
  Address = Address + FileName + dash + num_iter + Suffix;
  const char *ADDRESS = Address.c_str();

  int n1, n2;
  n1 = Vec2d.size();
  n2 = Vec2d[0].size();
  
  ofstream outfile;
  outfile.open(ADDRESS);
  
  for (int i=0; i<n1; i++) {
    for (int j=0; j<n2; j++) {
      outfile << setprecision(14) << Vec2d[i][j];
      if (j<n2-1)
	outfile<<"  ";
    }
    outfile << "\n";
  }  
  outfile.close();
  ADDRESS = new char; delete ADDRESS; ADDRESS = NULL;
}

void Output_1D_Vector(string FileName, vector<double> Vec)
{
  int VecSize = Vec.size();
  string Address = "./" + Output_Folder + "/";
  string Suffix = ".txt";
  Address = Address + FileName + Suffix;
  const char *ADDRESS = Address.c_str();

  ofstream outfile;
  outfile.open(ADDRESS);

  for (int i = 0; i <VecSize; i++)
    {
      outfile << setprecision(14)<< i <<"     "<< Vec[i] << endl;
    }
  outfile.close();
  ADDRESS = new char; delete ADDRESS; ADDRESS = NULL;
}

void Save_data_to_folder()
{  
  string execute1 = "exec rm ./OutputFiles/*";
  string execute2 = "cp -r "+Output_Folder+"/ OutputFiles/";
  const char *execline1 = execute1.c_str();
  const char *execline2 = execute2.c_str();
  system(execline1);
  system(execline2);

  execline1 = new char; delete execline1; execline1 = NULL;
  execline2 = new char; delete execline2; execline2 = NULL;
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
