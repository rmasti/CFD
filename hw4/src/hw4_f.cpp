/*
 * Computational Fluid Dynamics
 * Written By: Robert Masti
 * 1/27/2017
 * This file contains all function declarations that will be called in hw1_m.cpp which is the 
 * main cpp file. This is the function file, which has its function prototypes stored in the 
 * header file.
 */
#include "hw4.hpp" //structure templates and func prototypes



void compute_F_roe(
    // This function computes the F fluxes using Roe's flux diff scheme with Roe-averaged vars
    MatrixXd* F,                      //output - Flux at all interface of interest N+1     
    MatrixXd* V_L,                    //input - Left state primitive variables
    MatrixXd* V_R,                    //input - Right state primitive variables
    constants C)                      //input - constants C for gamma maybe
{
  // Because of the difficulty of handling expansion fan create small epsilon to asymptote to
  double eps = 0.1;

  // Declare roe averaged vars
  double rho_Roe; // rho avg
  double u_Roe; // u avg
  double v_Roe; // v avg
  double ht_Roe; // ht avg
  double a_Roe; // sound speed avg
  double R_Roe; // used to get rho_Roe

  // Declare LR ht
  double ht_L, ht_R;
  // Declare delta vars
  double drho, du, dv, dP;
  // define absolute value of lambda
  double abs_lambda1_Roe, abs_lambda2_Roe, abs_lambda3_Roe, abs_lambda4_Roe; 
  // define r's which go with the characteristic variables //4 vars long 
  double r1_Roe[neq], r2_Roe[neq], r3_Roe[neq], r4_Roe[neq];
  // define characteristic vars
  double dw1, dw2, dw3, dw4;

  // Write flux contributions L and R states
  double F1_L, F2_L, F3_L, F4_L;
  double F1_R, F2_R, F3_R, F4_R;


  // Loop over all of the interfaces of interest
  int row=0;
  for (int i = 0; i < F[frhouid].cols(); i++)
  {
    // get roe averaged roe
    R_Roe = sqrt(V_R[rhoid](row,i) / V_L[rhoid](row,i));
    rho_Roe = R_Roe*V_L[rhoid](row,i);

    // get roe averaged u
    u_Roe = (R_Roe*V_R[uid](row,i) + V_L[uid](row,i)) / (R_Roe + 1);

    // get roe averaged v
    // ...

    // get roe averaged ht
    ht_L = (C.gamma/(C.gamma - 1))*V_L[pid](row,i) / V_L[rhoid](row,i) 
      + 0.5*V_L[uid](row,i)*V_L[uid](row,i);
    ht_R = (C.gamma/(C.gamma - 1))*V_R[pid](row,i) / V_R[rhoid](row,i) 
      + 0.5*V_R[uid](row,i)*V_R[uid](row,i);
    ht_Roe = (R_Roe*ht_R + ht_L) / (R_Roe + 1);

    // get roe averaged soundspeed 
    a_Roe = sqrt( (C.gamma - 1)*(ht_Roe - 0.5*u_Roe*u_Roe));

    // compute eigen values
    abs_lambda1_Roe = abs(u_Roe);
    abs_lambda2_Roe = abs(u_Roe + a_Roe);
    abs_lambda3_Roe = 0; // v lambda uknown
    abs_lambda4_Roe = abs(u_Roe - a_Roe);

    // Apply expansion fan fix with epsilons
    if (abs_lambda1_Roe < 2*eps*a_Roe) 
      abs_lambda1_Roe = abs_lambda1_Roe*abs_lambda1_Roe/(4.0*eps*a_Roe) + eps*a_Roe;
    if (abs_lambda2_Roe < 2*eps*a_Roe) 
      abs_lambda2_Roe = abs_lambda2_Roe*abs_lambda2_Roe/(4.0*eps*a_Roe) + eps*a_Roe;
    if (abs_lambda3_Roe < 2*eps*a_Roe) 
      abs_lambda3_Roe = 0;
    if (abs_lambda4_Roe < 2*eps*a_Roe) 
      abs_lambda4_Roe = abs_lambda4_Roe*abs_lambda4_Roe/(4.0*eps*a_Roe) + eps*a_Roe;

    // get the right eigen vectors
    r1_Roe[0] = 1.0; //cont
    r1_Roe[1] = u_Roe; // xmtm
    r1_Roe[2] = 0;  // ymtm
    r1_Roe[3] = 0.5*u_Roe*u_Roe; // energy

    r2_Roe[0] = (0.5*rho_Roe/a_Roe)*(1.0);
    r2_Roe[1] = (0.5*rho_Roe/a_Roe)*(u_Roe + a_Roe);
    r2_Roe[2] = (0.5*rho_Roe/a_Roe)*(0.0);
    r2_Roe[3] = (0.5*rho_Roe/a_Roe)*(ht_Roe + u_Roe*a_Roe);
 
    r3_Roe[0] = 0.0;
    r3_Roe[1] = 0.0;
    r3_Roe[2] = 0.0;
    r3_Roe[3] = 0.0;

    r4_Roe[0] = -(0.5*rho_Roe/a_Roe)*(1.0);
    r4_Roe[1] = -(0.5*rho_Roe/a_Roe)*(u_Roe - a_Roe);
    r4_Roe[2] = -(0.5*rho_Roe/a_Roe)*(0.0);
    r4_Roe[3] = -(0.5*rho_Roe/a_Roe)*(ht_Roe - u_Roe*a_Roe);

    // get the variation
    drho = V_R[rhoid](row,i) - V_L[rhoid](row,i);
    du = V_R[uid](row,i) - V_L[uid](row,i);
    dv = V_R[vid](row,i) - V_L[vid](row,i);
    dP = V_R[pid](row,i) - V_L[pid](row,i);

    // find characteristic vars
    dw1 = drho - dP/(a_Roe*a_Roe);
    dw2 = du + dP/(rho_Roe*a_Roe);
    dw3 = 0.0;
    dw4 = du - dP/(rho_Roe*a_Roe);
    
    // calc first order fluxes left and right states
    // cont
    F1_L = V_L[rhoid](row,i)*V_L[uid](row,i);
    F1_R = V_R[rhoid](row,i)*V_R[uid](row,i);
    // x-mtm
    F2_L = V_L[rhoid](row,i)*V_L[uid](row,i)*V_L[uid](row,i) + V_L[pid](row,i);
    F2_R = V_R[rhoid](row,i)*V_R[uid](row,i)*V_R[uid](row,i) + V_R[pid](row,i);
    // y-mtm
    F3_L = 0.0;
    F3_R = 0.0;
    // energy
    F4_L = V_L[rhoid](row,i)*V_L[uid](row,i)*ht_L;
    F4_R = V_R[rhoid](row,i)*V_R[uid](row,i)*ht_R;

    // Compute fluxes
    // cont
    F[frhouid](row,i) = 0.5*(F1_L + F1_R)
      -0.5*(abs(abs_lambda1_Roe)*dw1*r1_Roe[0] 
        + abs(abs_lambda2_Roe)*dw2*r2_Roe[0]
        + abs(abs_lambda3_Roe)*dw3*r3_Roe[0]
        + abs(abs_lambda4_Roe)*dw4*r4_Roe[0]
        );

    // x-mtm
    F[frhouuid](row,i) = 0.5*(F2_L + F2_R) 
      -0.5*(abs(abs_lambda1_Roe)*dw1*r1_Roe[1] 
        + abs(abs_lambda2_Roe)*dw2*r2_Roe[1]
        + abs(abs_lambda3_Roe)*dw3*r3_Roe[1]
        + abs(abs_lambda4_Roe)*dw4*r4_Roe[1]
        );

    // y-mtm
    F[frhovid](row,i) = 0.5*(F3_L + F3_R) 
      -0.5*(abs(abs_lambda1_Roe)*dw1*r1_Roe[2] 
        + abs(abs_lambda2_Roe)*dw2*r2_Roe[2]
        + abs(abs_lambda3_Roe)*dw3*r3_Roe[2]
        + abs(abs_lambda4_Roe)*dw4*r4_Roe[2]
        );

    // x-mtm
    F[frhouhtid](row,i) = 0.5*(F4_L + F4_R)
      -0.5*(abs(abs_lambda1_Roe)*dw1*r1_Roe[3] 
        + abs(abs_lambda2_Roe)*dw2*r2_Roe[3]
        + abs(abs_lambda3_Roe)*dw3*r3_Roe[3]
        + abs(abs_lambda4_Roe)*dw4*r4_Roe[3]
        );
  }
}

void compute_F_vanleer(
    // This function will use the right and left state primitive variables to determine fluxes
    MatrixXd* F,                      //output - Flux at all interfaces of interest N+1        
    MatrixXd* V_L,                    //input - Left state at interfaces of interest
    MatrixXd* V_R,                    //input - Right state at interfaces of interest
    constants C)                      //input - constants C
{
  double M_L, M_R, M_pos, M_neg; // Mach left right then pos neg
  double beta_L, beta_R; // beta left right using in final form
  double alpha_pos, alpha_neg; // alpha pos neg used in calc c
  double c_pos, c_neg; // c is used in flux addition convective
  double a_L, a_R; // sound speed
  double D_pos, D_neg; // pressure flux constants
  double Pbar_pos, Pbar_neg; // pressure bar for mtm?
  double ht_L, ht_R; // enthalpy total left and right states
  // Define Flux Contributions
  double Fc1, Fc2, Fc3, Fc4; // frhouid, frhouuid, frhovid, frhouhtid //convective flux
  double Fp1, Fp2, Fp3, Fp4; // frhouid, frhouuid, frhovid, frhouhtid //Pressure flux

  int row = 0;

  // loop over and fill in flux
  for (int i = 0; i < F[frhouid].cols(); i++)
  {
    // Find the soundspeeds
    a_L = compute_soundspeed(C.gamma, V_L[pid](row,i),V_L[rhoid](row,i));
    a_R = compute_soundspeed(C.gamma, V_R[pid](row,i),V_R[rhoid](row,i));

    // Find the Mach numbers left and right state
    M_L = V_L[uid](row,i)/a_L;
    M_R = V_R[uid](row,i)/a_R;

    // Find the pos and negative mach
    M_pos = 0.25*(M_L + 1)*(M_L + 1);
    M_neg = -0.25*(M_R - 1)*(M_R - 1);

    // Find Beta left and right state
    beta_L = -mymax(0,1 - int(abs(M_L)));
    beta_R = -mymax(0,1 - int(abs(M_R)));

    // Find alpha pos and neg
    alpha_pos = 0.5*(1 + SIGN(1,M_L));
    alpha_neg = 0.5*(1 - SIGN(1,M_R));

    // Find c pos neg for use in flux finding
    c_pos = alpha_pos*(1 + beta_L)*M_L - beta_L*M_pos;
    c_neg = alpha_neg*(1 + beta_R)*M_R - beta_R*M_neg;

    // Find the Pbar pos and neg
    Pbar_pos = M_pos*(-M_L + 2);
    Pbar_neg = M_neg*(-M_R - 2);

    // Find the D pos and neg
    D_pos = alpha_pos*(1 + beta_L) - beta_L*Pbar_pos;
    D_neg = alpha_neg*(1 + beta_R) - beta_R*Pbar_neg;

    // Find ht L and R states
    ht_L = (C.gamma/(C.gamma - 1))*V_L[pid](row,i) / V_L[rhoid](row,i) 
      + 0.5*V_L[uid](row,i)*V_L[uid](row,i);
    ht_R = (C.gamma/(C.gamma - 1))*V_R[pid](row,i) / V_R[rhoid](row,i) 
      + 0.5*V_R[uid](row,i)*V_R[uid](row,i);

    // Find Convective Flux Contributions
    // cont
    Fc1 = V_L[rhoid](row,i)*a_L*c_pos + V_R[rhoid](row,i)*a_R*c_neg;
    // x-mtm
    Fc2 = V_L[rhoid](row,i)*a_L*c_pos*V_L[uid](row,i)+V_R[rhoid](row,i)*a_R*c_neg*V_R[uid](row,i);
    // y-mtm
    Fc3 = 0.0;
    // ht 
    Fc4 = V_L[rhoid](row,i)*a_L*c_pos*ht_L + V_R[rhoid](row,i)*a_R*c_neg*ht_R;

    // Find Pressure Flux Contributions
    // cont
    Fp1 = 0.0;
    // x-mtm
    Fp2 = D_pos*V_L[pid](row,i) + D_neg*V_R[pid](row,i);
    // y-mtm
    Fp3 = 0.0;
    // ht
    Fp4 = 0.0;

    // Compute the flux at the interfaces
    F[frhouid](row,i) = Fc1 + Fp1;
    F[frhouuid](row,i) = Fc2 + Fp2;
    F[frhovid](row,i) = Fc3 + Fp3;
    F[frhouhtid](row,i) = Fc4 + Fp4;
  }
}


void compute_upwind_VLR(
    // Compute the upwind Left and right primvar states 
    MatrixXd* V_L,                    //output - Left primvar state at interfaces of interest           
    MatrixXd* V_R,                    //output - Right primvar state at interfaces of interest
    MatrixXd* V,                      //input - Primvar at all cell centers
    constants C)                     //input - Choose limiter                     
{
  // Create the psi_positive and psi_negative and define at all interfaces So N+ghost+1
  //
  int row = 0;
  MatrixXd*  Psi_Pos=new MatrixXd[neq];
  MatrixXd*  Psi_Neg=new MatrixXd[neq];
  for (int eq = 0 ; eq < neq; eq++)
  {
    Psi_Pos[eq].resize(1,V[rhoid].cols()+1); //Every interface can have a psi value*****
    Psi_Neg[eq].resize(1,V[rhoid].cols()+1); //Also can have negative values
  }
  // At this point Psi's have the correct size, but will just need to be filled in now
  compute_psi_pn(Psi_Pos, Psi_Neg, V, C);

  // Now the psi values are known at all interface values now fill in LR states
  for (int i = 0 ; i < V_L[rhoid].cols(); i++)
  {
    int i_cells = (num_ghost_cells-1) + i; // start at i = 2 where i+1/2 starts 
    int i_int = num_ghost_cells+i; // starts at the 3rd interface
    for (int eq = 0 ; eq < neq; eq++)
    {    
      // Note i_cells + 1 is the index for the i+1/2 value of psi
      // Left state at the interface
      V_L[eq](row,i) = V[eq](row,i_cells) + 0.25*upwind_epsilon*(
          (1.0-upwind_kappa)*Psi_Pos[eq](row,i_int-1)*(V[eq](row,i_cells)-V[eq](row,i_cells-1))
          + (1.0+upwind_kappa)*Psi_Neg[eq](row,i_int)*(V[eq](row,i_cells+1)-V[eq](row,i_cells))
          );
      // Right state at the interface
      V_R[eq](row,i) = V[eq](row,i_cells+1) - 0.25*upwind_epsilon*(
          (1.0-upwind_kappa)*Psi_Neg[eq](row,i_int+1)*(V[eq](row,i_cells+2)-V[eq](row,i_cells+1))
          + (1.0+upwind_kappa)*Psi_Pos[eq](row,i_int)*(V[eq](row,i_cells+1)-V[eq](row,i_cells))
          );
    }
  }
}


void compute_psi_pn(
    // Fill in the Pos and Neg Psi at all interface values possible 2 - last-2
    // This also applies the limiters to the psi values
    MatrixXd* Psi_Pos,                //ouput - Psi_Pos at interfaces                   
    MatrixXd* Psi_Neg,                //output - Psi_Neg at interfaces
    MatrixXd* V,                      //input - Used to find slopes
    constants C)                      //input - Choose the limiter
{

  int row = 0;
  double delta = 1.0e-6; // Small number to avoid division by 0
  // Loop over the cell centroid indexes then just add 1 for the i+1/2 interface index
  for (int i = num_ghost_cells-2 ; i <= N+num_ghost_cells; i++)
  {
    // Loop over all of the primvars
    for (int eq = 0; eq < neq; eq++)
    {
      // Calculate the denom which is the same for pos and neg
      double denom = V[eq](row,i+1) - V[eq](row,i);
      // Ensure not division by 0
      denom = SIGN( mymax(abs(denom), delta), denom);
      // Calculate the slopes (r) for both pos and neg dir
      double r_p = (V[eq](row,i+2) - V[eq](row, i+1))/denom;
      double r_n = (V[eq](row,i) - V[eq](row, i-1))/denom;
      // Apply Limiters
      if ( C.limiter == 0 )  //no limiters
      {
        Psi_Pos[eq](row,i+1) = 1.0;
        Psi_Neg[eq](row,i+1) = 1.0;
      }
      if ( C.limiter == 1 )  //Van leer
      {
        Psi_Pos[eq](row,i+1) = (r_p + abs(r_p)) / (1 + r_p);
        Psi_Neg[eq](row,i+1) = (r_n + abs(r_n)) / (1 + r_n);
      }
      if ( C.limiter == 2 ) //Van Albada 
      {
        // if negative reduce to first order
        if( r_p < 0)
          Psi_Pos[eq](row,i+1) = 0.0;
        else
          Psi_Pos[eq](row,i+1) = (r_p + r_p*r_p) / (1 + r_p*r_p);
        if( r_n < 0)
          Psi_Neg[eq](row,i+1) = 0.0;
        else
          Psi_Neg[eq](row,i+1) = (r_n + r_n*r_n) / (1 + r_n*r_n);
      }
    }
  }
}

void write_solution(
    //This function outputs all required variables into a 
    //single matrix for the whole 1d solution
    FILE*& file,                      //output - Values to file
    MatrixXd& xc,                     //input - xcoord at centers
    MatrixXd& Ac,                     //input - Area at centers
    MatrixXd* V,                      //input - Primitive vars at centers
    MatrixXd* U,                      //input - conservative vars at centers
    constants C)                      //input - mach number only
{
  int row = 0;
  //Loop over domain of interest only
  for(int i = num_ghost_cells; i < N+num_ghost_cells; i++)
  {
    //compute mach number
    double a = compute_soundspeed(C.gamma, V[pid](row,i),V[rhoid](row,i));
    double M = V[uid](row,i)/a; 
    //output the solution
    fprintf(file, "%e %e %e %e %e %e %e %e %e %e %e\n",xc(row,i), Ac(row,i), 
        V[rhoid](row,i), V[uid](row,i), V[vid](row,i), V[pid](row,i), M, 
        U[rhoid](row,i), U[rhouid](row,i), U[rhovid](row,i), U[rhoetid](row,i));
  }
}

void compute_norms(
    //This function takes in the residuals and finds the l2 norm
    MatrixXd& L2norm,                 //output - L2norm of all neq
    MatrixXd* Res)                    //input - Residuals on the whole domain
{
  int row = 0;
  for(int eq=0 ; eq < neq; eq++)
  {
    MatrixXd temp = Res[eq].cwiseProduct(Res[eq]);
    L2norm(row,eq) = sqrt(temp.sum()/N);
  }
}

void constoprim(
    //This function takes Uarray and converts it to a primitive variable
    MatrixXd* V,                      //output - Primitive var
    MatrixXd* U,                      //input - Conserved var
    constants C)                      //input - Constants
{
  V[rhoid] = U[rhoid];
  V[uid] = U[rhouid].cwiseProduct(U[rhoid].cwiseInverse()); // "rhou/rho"
  V[vid] = U[rhovid].cwiseProduct(U[rhoid].cwiseInverse()); // "rhov/rho"
  V[pid] = (C.gamma - 1.0)*(U[rhoetid] - 
      0.5*V[rhoid].cwiseProduct(V[uid].cwiseProduct(V[uid])) -
      0.5*V[rhoid].cwiseProduct(V[vid].cwiseProduct(V[vid]))); // "rhov/rho"
}
void iteration(
    //takes a time step iteration over the interior domain only!!!!!
    MatrixXd* U,                       //output - Updates Conserved Variable
    double& timestep,                  //output - timestep for keeping track of t
    MatrixXd* Res,                     //input - Residual of array size N
    MatrixXd& Ac,                      //input - Area at cell centers for 1d
    MatrixXd& Lambda_mcenter,          //input - Lambda max at center values for time step
    constants C)                       //input - constants for cfl
{

  MatrixXd vol = Ac*(dx); //vol now for N+2*numghost
  MatrixXd dtvec = C.cfl*dx*Lambda_mcenter.cwiseInverse();//get the dt vector from lambda's
  double dt = dtvec.minCoeff();//grab the minimum value of the matrix 
  timestep = dt;
  int row = 0;
  for(int i = 0 ; i < Res[rhoid].cols(); i++)
  {
    int i_cells = num_ghost_cells+i;//i_cells 3->13 
    for(int eq=0; eq < neq; eq++)
    {
      U[eq](row, i_cells) = U[eq](row,i_cells) - (dt/vol(row,i_cells))*Res[eq](row,i);
    }
  }
}

void compute_residual(
    //This function computes the residual at all N cells inside the domain
    MatrixXd* Res,                      //output - Residual array of size N
    MatrixXd* S,                        //input - Source array
    MatrixXd* F,                        //input - F flux size N+1
    MatrixXd& Ai)                       //input - Area at interfaces size N+2*numghost+1
{
  int row = 0;
  for(int i = 0; i < Res[rhoid].cols(); i++)
  {
    int i_int = num_ghost_cells + i; //for area vector
    for(int eq=0; eq < neq; eq++)
    {
      Res[eq](row,i) = (F[eq](row,i+1))*Ai(row,i_int+1) -  
        (F[eq](row,i))*Ai(row,i_int) - dx*S[eq](row,i);
    }
  }
}

void compute_source(
    // This function computes the source term for the rhou component of U in 1D
    MatrixXd* S,                        //output - Source Term calc at centers
    MatrixXd* V,                        //input - Prim var for pressures
    MatrixXd& xc)                       //input - xcenters for dAdx 
{
  int row = 0;
  for(int i=0; i < S[rhoid].cols(); i++)
  {
    int i_cells = i+num_ghost_cells;
    S[rhoid](row,i) = 0;
    S[rhouid](row,i) = V[pid](row,i_cells)*dAdx(xc(row,i_cells));//Pressure*dAdx
    S[rhovid](row,i) = 0;
    S[rhoetid](row,i) = 0;
  }
}
void compute_nu(
    // This function creates the matrix of nu's needed to compute epsilon
    MatrixXd& nu,                       //output - nu for use in epsilon NOTE N+4 size
    MatrixXd* V)                        //input - V for the pressures
{
  int row=0;
  nu(row,0) = 0.0;
  nu(row,nu.cols()-1) = 0.0;
  for(int i = 1; i < nu.cols()-1; i++)
  {
    nu(row,i) = abs((V[pid](row,i+1) - 2.0*V[pid](row,i) + V[pid](row,i-1)))/
      (V[pid](row,i+1) + 2.0*V[pid](row,i) + V[pid](row,i-1));
    //cout << nu(row,i)<< endl;
  }
}
void add_artificial_viscosity(
    // This function will compute dflux using jameson damping for all interfaces of interest
    MatrixXd* F,                         //output - F fluxes at all interfaces
    MatrixXd* V,                         //input - Prim var for P's
    MatrixXd* U,                         //input - Cons var for D1 and D3
    MatrixXd& Lambda_minterface)         //input - Lambda max at interfaces of interest
{

  int row = 0;
  // Fill the nu values only needed for art visc
  // NOTE that it is the same size as the V vectors
  // because they can use the same index in the loop below
  // nu is filled from 2 -> N+2*numghost-2
  MatrixXd nu(1,N+2*num_ghost_cells);
  compute_nu(nu, V);
  // Loop over the lambda of interest values
  for(int i = 0; i < Lambda_minterface.cols(); i++)
  {
    int i_cells = (num_ghost_cells-1) + i; //i_cells starts at 2 - > N+ghost 
    // Find the maximum nu->epsion2, and epsion 4
    double epsilon2 = kappa2*mymax(
        mymax(nu(row,i_cells-1),nu(row,i_cells)),
        mymax(nu(row,i_cells+1),nu(row,i_cells+2)));

    double epsilon4 = mymax(0, (kappa4 - epsilon2));
    // write lam for brevity in the long calcs
    double lam = Lambda_minterface(row,i);

    // Compute D1 D3 and d
    MatrixXd D1(1,neq);
    MatrixXd D3(1,neq);
    for(int eq = 0; eq < neq ; eq++)
    {
      D1(row,eq) = lam*epsilon2*(U[eq](row,i_cells+1) - U[eq](row,i_cells));
      D3(row,eq) = lam*epsilon4*(U[eq](row,i_cells+2) - 3.0*U[eq](row,i_cells+1)+
          3.0*U[eq](row,i_cells) - U[eq](row,i_cells-1));
      //combine D3 and D1
      // add d flux to the interface flux as is so it takes the d contributions and adds it to F
      F[eq](row,i) = F[eq](row,i) + D3(row,eq) - D1(row,eq); 
    }
  }
}
void compute_F_jameson(
    //This function takes in the U values at all the interfaces of interest and converts them to
    //the flux F
    MatrixXd* F,                                    //output - Flux from Uint
    MatrixXd* Uinterface,                           //input - U at i+1/2 for domain of Interest
    constants C)                                    //input - constants for conversion
{
  int row = 0;
  for (int i=0; i < Uinterface[rhoid].cols(); i++)
  {
    double rho = Uinterface[rhoid](row,i);
    double rhou = Uinterface[rhouid](row,i);
    double rhov = Uinterface[rhovid](row,i);
    double rhoet = Uinterface[rhoetid](row,i);

    //compute first element rho*u
    F[frhouid](row,i) = rhou; //first element is rhou

    //compute second element rho*u*u+p
    F[frhouuid](row,i) = 0.5*(3.0-C.gamma)*(rhou*rhou/rho) + 
      (C.gamma-1.0)*rhoet;//second element is rhouu + p

    //compute third element rho*v for F NOT G
    F[frhovid](row,i) = rhov;

    //compute fourth element rho*u*ht
    F[frhouhtid](row,i) = (rhoet*rhou/rho) + (rhou/rho)*
      (C.gamma - 1.0)*(rhoet - 0.5*rhou*rhou/rho);//fourth element is rhouht
  }
}

void reconstruct(
    //this function recunstructs cell values to the interfaces ONLY OF INTEREST N+1 no ghost
    MatrixXd* Uinterface,                           //output - U at interfaces of interest
    MatrixXd& Lambda_minterface,                    //output - lambda at interfaces of interest
    MatrixXd* U,                                    //input - U at cell values
    MatrixXd& Lambda_mcenter)                       //input - lambda at interfaces
{
  int row = 0;
  for(int i = 0; i < Lambda_minterface.cols(); i++)
  {
    int i_cells = num_ghost_cells + i; //from 3 -> N+3 if num ghost = 3
    Lambda_minterface(row,i) = 0.5*(Lambda_mcenter(i_cells) + Lambda_mcenter(i_cells-1));
    // reconstrunct conservative variables
    for(int eq = 0; eq < neq; eq++)
    {
      Uinterface[eq](row,i) = 0.5*(U[eq](row,i_cells) + U[eq](row,i_cells-1));
    }
  }
}

void compute_lambda(
    //This function computes lambdamax at centers by avg lambda from centers
    MatrixXd& Lambda_mcenter,                      //output - Lambda at centers
    MatrixXd* V,                                   //input - Prim var at centers
    constants C)                                   //input - constants
{
  int row = 0;
  for(int i = 0; i < Lambda_mcenter.cols(); i++)
  {
    Lambda_mcenter(row,i) = abs(V[uid](row,i)) + 
      compute_soundspeed(C.gamma, V[pid](row,i), V[rhoid](row,i));
  }
}

void primtocons(
    //this function converts primitive variables to conservative variables
    MatrixXd* U,                                   //output - Conservative Variables
    MatrixXd* V,                                   //input - Primitive Variables
    constants C)                                   //input - Constants
{
  U[rhoid] = V[rhoid];//U1 = V1
  U[rhouid] = V[rhoid].cwiseProduct(V[uid]); //U2 = rhou = V1*V2
  U[rhovid] = V[rhoid].cwiseProduct(V[vid]); //U3 = rhov = V1*V3
  //Compute energy U4 = V4/(gamma-1.0) + 0.5*V1*V2*V2 + 0.5*V1*V3*V3
  U[rhoetid] = V[pid]/(C.gamma - 1.0) 
    + 0.5*V[rhoid].cwiseProduct(V[uid].cwiseProduct(V[uid]));
}

void outflow_boundary_condition(
    //this specifies outflow only for subsonic outlet
    MatrixXd* V,                                   //output - primitive variables
    constants C)                                   //input - constants
{
  if(C.outflow == false) //then subsonic and requires extra treatment
  {
    int row = 0;
    int end = V[pid].cols()-1;//the very last cell
    double P_middle = C.pb; //edge boundary or -> prev cell
    double P_left = V[pid](row, end-num_ghost_cells); //last interior cell
    for(int i = N+num_ghost_cells; i < V[pid].cols() ;i++)
    {
      V[pid](row, i) = C.pb;
      P_left = V[pid](row, i-1);
      P_middle = V[pid](row, i); 
    }
  }
}

void inflow_boundary_condition(
    //this specifies inflow only
    MatrixXd* V,                                    //output - primitive variables
    constants C)                                    //input - constants
{
  int row = 0;  
  double M_1 = V[uid](row,num_ghost_cells)/compute_soundspeed(C.gamma, V[pid](row,num_ghost_cells), V[rhoid](row,num_ghost_cells)); //u/a
  double M_2 = V[uid](row,num_ghost_cells+1)/compute_soundspeed(C.gamma, V[pid](row,num_ghost_cells+1), V[rhoid](row,num_ghost_cells+1)); //u/a
  double M_0 = 0;
  //NEED TO extrapolate mach number from the interior out to the ghost cells
  for(int i=0; i< num_ghost_cells; i++)
  {
    M_0 = 2.0*M_1 - M_2;
    if (M_0 < 0.001){
      M_0 = 0.001;
      cout << "Positivity Preserving Activated!!!" << endl;
    }
    //now with the mach number apply Mto prims with the P0 and T0 specified
    int ib = (num_ghost_cells - 1) - i; 
    Mtoprim(V[rhoid](row,ib), V[uid](row,ib), V[vid](row,ib), V[pid](row,ib), M_0, C);
    M_2 = M_1;
    M_1 = M_0;
  }
}

void set_boundary_conditions(
    //This function will apply boundary conditions to the primitive var array
    MatrixXd* V,                                      //output - primitive variables
    constants C)                                      //input - constants
{
  inflow_boundary_condition(V, C);
  outflow_boundary_condition(V,C);
}

void extrapolate_to_ghost(MatrixXd* Varr)
  //This function takes and array of matrix an extrapolates to the ghost cells
{
  int row = 0;//for 2d
  int end = Varr[rhoid].cols()-1;//helps with determining indexes
  for(int i = 0; i < num_ghost_cells; i++)
  {
    int ileft = num_ghost_cells - i; //Start at the 1st interior cell and extrap to ghost
    int iright = (end - num_ghost_cells) + i;//Find the last interior cell
    for(int eq=0; eq < neq; eq++)
    {
      // Left extrap
      Varr[eq](row, ileft-1) = 2.0*Varr[eq](row, ileft) - Varr[eq](row, ileft+1);//rho
      // Right extrap
      Varr[eq](row, iright+1) = 2.0*Varr[eq](row, iright) - Varr[eq](row, iright-1);//rho
    }
  }
}

void initialize(
    //This function initialize the primitive variable at all the cell centers throughout 
    //the domain %%%%%%% Including Ghost Cells %%%%%%%%%%
    MatrixXd* V,                            //output - The primitive variables at centers 
    MatrixXd &Mc,                           //input - mach number at cell centers
    constants &C)                           //input - Constants is passed for conversions
{
  //Define values used as inputs for the Mtoprim function
  for (int row=0; row < Mc.rows(); row++) 
  {
    for(int col = 0; col < Mc.cols(); col++)
    {
      Mtoprim(V[rhoid](row,col), V[uid](row,col),
          V[vid](row,col), V[pid](row,col), Mc(row,col), C);// fill in V1-4
    }
  }
}

void Mtoprim(
    //During Initialization the Mach number is initialized first which then 
    //needs to be converted to a primitive variable that is what this func does
    double& V1,                                //output - Rho
    double& V2,                                //output - u-velocity
    double& V3,                                //output - v-velocity
    double& V4,                                //output - Pressure 
    double& M,                                 //input - Mach number 
    constants& C)                              //input - Constants
{
  double psi = 1.0 + 0.5*(C.gamma - 1.0) * M*M;//unitless
  double T = C.T0/psi;//Kelvin Temperature calc 
  V4 = (C.p0)/pow(psi,C.gamma/(C.gamma - 1.0));//Pa Pressure calc
  V1 = V4/(R*T); //kg/m^3 Rho calc
  double a = compute_soundspeed(C.gamma, V4, V1);//compute sound speed
  V2 = a*M; //m/s u calc
  V3 = 0.0; //m/s v calc 0.0 because in 1D
}

void set_geometry(
    //This function sets the x interface coordinates, and also the x center coordinates
    MatrixXd& xcenter,                        //output - xcoord cell center            
    MatrixXd& xinterface)                     //output - xcoord cell interfaces
{
  double x_int_left = xmin-num_ghost_cells*dx;// left most x coord
  for (int i=0; i < xcenter.cols() ; i++)//.cols() returns the number of columns 
  {
    // create the coordinate
    double x_i = x_int_left + dx*i;
    double x_c = x_i + dx/2.0;
    // write to the matrix
    xinterface(0,i) = x_i;
    xcenter(0,i) = x_c;
  }
  // Add the last interface value
  xinterface(0,xinterface.cols()-1) = xmax+num_ghost_cells*dx;
}

double A_x(double x)
  //This function takes in x-values and returns the area
{
  if (x < xmin_dom || x > xmax_dom) return 1.0;
  else return 0.2 + 0.4*(1.0 + sin(M_PI*(x-0.5)));//function that governs Area for nozzle
}

double M_xinitial(double x)
  //This function initializes the mach number everywhere in the domain
{
  double slope = (1.9 - 0.1)/(xmax_dom - xmin_dom);
  double b = 0.1-slope*xmin_dom;
  double out; 
  if (x < xmin_dom) out = slope*xmin_dom+b;
  else if (x > xmax_dom ) out = slope*xmax_dom+b;
  else out = slope*x+b;
  return out;
}

double compute_soundspeed(
    //This function computes the sound speed for p, gamma, and rho
    double gamma, 
    double p, 
    double rho) 
{
  return sqrt(gamma*p/rho);
}

double dAdx(double x)
{
  if (x < xmin_dom || x > xmax_dom) return 0.0;
  else return 0.4*M_PI*cos(M_PI*(x-0.5));//derivative of func above
}

double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}


void exactsol(
    //This function computes the exact solution
    MatrixXd* V,                               //output - Primitive variables including ghost
    MatrixXd& Ac,                              //input - Area at centers
    constants C)                               //input - Constants  for conversion
{
  double A_t = C.A_t;//m^2
  double tol = C.tol;//unitless
  double phi;//unitless
  double dFdM;//unitless
  double Mnew;//unitless
  double M;

  int row = 0;
  for (int i = 0; i < V[rhoid].cols(); i++)
  {
    double A_bar = Ac(i)/A_t;//unitless
    double F = 2.0;//unitless
    double dM=100000.0;//unitless
    if (i<((N+2*num_ghost_cells)/2)) M=0.5;
    else M=5.0;
    while (fabs(dM) > tol)
    {
      phi = (2.0/(C.gamma + 1.0)) * (1.0 + ((C.gamma - 1.0)/2.0) * M*M);//unitless
      F = pow(phi,(C.gamma + 1.0)/(C.gamma - 1.0)) - A_bar*A_bar * M*M;//unitless
      dFdM = 2.0 * M * (pow(phi,2.0/(C.gamma - 1.0)) - A_bar*A_bar);//unitless
      Mnew = M - 0.1*(F/dFdM); //unitless 
      dM = (Mnew-M);//Used to check convergence
      M = Mnew;
    }
    Mtoprim(V[rhoid](row,i), V[uid](row,i), V[vid](row,i), V[pid](row,i), M, C);
  }
  //return the primary variables V once the Mach number is known.
}
