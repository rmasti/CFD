/*
 * Computational Fluid Dynamics
 * Written ByRobert Masti
 * 1/27/2017
 * This file contains all function declarations that will be called in fp_m.cpp which is the 
 * main cpp file. This is the function file, which has its function prototypes stored in the 
 * header file.
 */
#include "fp.hpp" //structure templates and func prototypes and libs

void computeFluxVL(
    // This function will compute the van leer flux for arb left and right
    // it will then update flux 
    MatrixXd* Flux,        // output - flux     
    MatrixXd nxhat,        // input - norm vec xhat
    MatrixXd nyhat,        // input - norm vec yhat
    MatrixXd* V_Left,         // input - left states
    MatrixXd* V_Right,         // input - right states
    constants C            // input - constants for MHD
    )
{
  // Grab the matrices for the sake of cleanliness
  MatrixXd rho_L = V_Left[rhoid];
  MatrixXd rho_R = V_Right[rhoid];

  MatrixXd u_L = V_Left[uid];
  MatrixXd u_R = V_Right[uid];

  MatrixXd v_L = V_Left[vid];
  MatrixXd v_R = V_Right[vid];

  MatrixXd p_L = V_Left[pid];
  MatrixXd p_R = V_Right[pid];

  // define constants
  double M_L, M_R, M_p, M_n;
  double beta_L, beta_R;
  double alpha_p, alpha_n;
  double c_p, c_n;
  double a_L, a_R;
  double U_L, U_R;
  double D_p, D_n;
  double p_bar_p, p_bar_n;
  double ht_L, ht_R; 
  double Fc[NEQ], Fp[NEQ];

  int ni = Flux[0].cols();
  int nj = Flux[0].rows();

  if (nxhat.rows() != nj || nxhat.cols() != ni)
  {
    cerr << "ERROR: Dimension Mismatch in VL Flux?!?!" << endl;
    exit(1);
  } 
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      a_L = sqrt(GAMMA*p_L(j,i) / rho_L(j,i)); 
      a_R = sqrt(GAMMA*p_R(j,i) / rho_R(j,i)); 

      U_L = u_L(j,i)*nxhat(j,i) + v_L(j,i)*nyhat(j,i);
      U_R = u_R(j,i)*nxhat(j,i) + v_R(j,i)*nyhat(j,i);

      M_L = U_L/a_L;
      M_R = U_R/a_R;

      M_p = 0.25*(M_L + 1)*(M_L + 1);
      M_n = -0.25*(M_R - 1)*(M_R - 1);

      beta_L = -mymax(0, 1 - int(abs(M_L)));
      beta_R = -mymax(0, 1 - int(abs(M_R)));

      alpha_p = 0.5*(1 + SIGN(1, M_L));
      alpha_n = 0.5*(1 - SIGN(1, M_R));

      c_p = alpha_p*(1+beta_L)*M_L - beta_L*M_p;
      c_n = alpha_n*(1+beta_R)*M_R - beta_R*M_n;

      p_bar_p = M_p*(-M_L + 2);
      p_bar_n = M_n*(-M_R - 2);

      D_p = alpha_p*(1+beta_L) - beta_L*p_bar_p;
      D_n = alpha_n*(1+beta_R) - beta_R*p_bar_n;

      ht_L = (GAMMA/(GAMMA-1))*p_L(j,i)/rho_L(j,i) 
        + 0.5*( u_L(j,i)*u_L(j,i) + v_L(j,i)*v_L(j,i));
      ht_R = (GAMMA/(GAMMA-1))*p_R(j,i)/rho_R(j,i) 
        + 0.5*( u_R(j,i)*u_R(j,i) + v_R(j,i)*v_R(j,i));

      Fc[0] = rho_L(j,i)*a_L*c_p + rho_R(j,i)*a_R*c_n;
      Fc[1] = rho_L(j,i)*a_L*c_p*u_L(j,i) + rho_R(j,i)*a_R*c_n*u_R(j,i);
      Fc[2] = rho_L(j,i)*a_L*c_p*v_L(j,i) + rho_R(j,i)*a_R*c_n*v_R(j,i);
      Fc[3] = rho_L(j,i)*a_L*c_p*ht_L + rho_R(j,i)*a_R*c_n*ht_R;

      Fp[0] = 0.0;
      Fp[1] = D_p*nxhat(j,i)*p_L(j,i) + D_n*nxhat(j,i)*p_R(j,i);
      Fp[2] = D_p*nyhat(j,i)*p_L(j,i) + D_n*nyhat(j,i)*p_R(j,i);
      Fp[3] = 0.0;

      for (int eq = 0; eq < NEQ; eq++)
        Flux[eq](j,i) = Fc[eq] + Fp[eq];
    }
  }
}

void computeFluxRoe(
    // This function will compute the flux for an arbitrary 
    // Left and right state direction which means it will work
    // for top and bottom, it will modify F1-F4,
    MatrixXd* Flux,        // output - flux component cont 
    MatrixXd nxhat,        // input - normal vector xhat comp         
    MatrixXd nyhat,        // input - normal vector yhat comp
    MatrixXd* V_Left,         // input - Left state
    MatrixXd* V_Right,         // input - Right state
    constants C            // input - constants for MHD maybe
    )
{

  // Grab the matrices for the sake of cleanliness
  MatrixXd rho_L = V_Left[rhoid];
  MatrixXd rho_R = V_Right[rhoid];

  MatrixXd u_L = V_Left[uid];
  MatrixXd u_R = V_Right[uid];

  MatrixXd v_L = V_Left[vid];
  MatrixXd v_R = V_Right[vid];

  MatrixXd p_L = V_Left[pid];
  MatrixXd p_R = V_Right[pid];

  // define roe vars
  double rho_Roe;
  double u_Roe;
  double v_Roe;
  double U_hat_Roe_R, U_hat_Roe_L;
  double U_hat_Roe;
  double ht_L, ht_R;
  double ht_Roe;
  double a_Roe;
  double R_Roe;

  double drho, du, dv, dp;
  double dw1, dw2, dw3, dw4;

  // eigen vals
  double lambda1, lambda2, lambda3, lambda4;
  double eps = 0.1; // correction term
  double r1[NEQ], r2[NEQ], r3[NEQ], r4[NEQ];

  double F_L[NEQ];
  double F_R[NEQ];

  double sum2ndOrder;

  int ni = u_L.cols();
  int nj = u_L.rows();

  if (nxhat.rows() != nj || nxhat.cols() != ni)
  {
    cerr << "ERROR: Dimension Mismatch in Roe Flux?!?!" << endl;
    exit(1);
  } 

  for (int j = 0; j < nj; j++)
  {
    for (int i = 0; i < ni; i++)
    {
      R_Roe = sqrt(rho_R(j,i)/rho_L(j,i));
      rho_Roe = R_Roe*rho_L(j,i);
      u_Roe = (R_Roe*u_R(j,i) + u_L(j,i)) / (R_Roe + 1);
      v_Roe = (R_Roe*v_R(j,i) + v_L(j,i)) / (R_Roe + 1);
      U_hat_Roe = u_Roe*nxhat(j,i) + v_Roe*nyhat(j,i);// get the speed

      ht_L = (GAMMA/(GAMMA-1))*p_L(j,i)/rho_L(j,i) 
        + 0.5*( u_L(j,i)*u_L(j,i) + v_L(j,i)*v_L(j,i));
      ht_R = (GAMMA/(GAMMA-1))*p_R(j,i)/rho_R(j,i) 
        + 0.5*( u_R(j,i)*u_R(j,i) + v_R(j,i)*v_R(j,i));
      ht_Roe = (R_Roe*ht_R + ht_L) / (R_Roe + 1);

      a_Roe = sqrt((GAMMA-1)*( ht_Roe - 0.5*(u_Roe*u_Roe + v_Roe*v_Roe)));

      lambda1 = abs(U_hat_Roe); // speed 
      lambda2 = abs(U_hat_Roe); // speed
      lambda3 = abs(U_hat_Roe + a_Roe); // u+a 1D
      lambda4 = abs(U_hat_Roe - a_Roe); // u-a 1D

      if (lambda1 < 2*eps*a_Roe)
      {
        lambda1 = lambda1*lambda1/(4*eps*a_Roe) + eps*a_Roe;
        lambda2 = lambda2*lambda2/(4*eps*a_Roe) + eps*a_Roe;
      }
      if (lambda3 < 2*eps*a_Roe)
        lambda3 = lambda3*lambda3/(4*eps*a_Roe) + eps*a_Roe;
      if (lambda4 < 2*eps*a_Roe)
        lambda4 = lambda4*lambda4/(4*eps*a_Roe) + eps*a_Roe;

      // fill in the right eigen vectors
      r1[0] = 1.0;
      r1[1] = u_Roe;
      r1[2] = v_Roe;
      r1[3] = 0.5*(u_Roe*u_Roe + v_Roe*v_Roe);

      r2[0] = 0.0;
      r2[1] = rho_Roe * nyhat(j,i);
      r2[2] = -rho_Roe * nxhat(j,i);
      r2[3] = rho_Roe*(u_Roe*nyhat(j,i) - v_Roe*nxhat(j,i));

      r3[0] = (0.5*rho_Roe/a_Roe);
      r3[1] = (0.5*rho_Roe/a_Roe) * ( u_Roe + a_Roe*nxhat(j,i) );
      r3[2] = (0.5*rho_Roe/a_Roe) * ( v_Roe + a_Roe*nyhat(j,i) );
      r3[3] = (0.5*rho_Roe/a_Roe) * ( ht_Roe + a_Roe*U_hat_Roe );

      r4[0] = (-0.5*rho_Roe/a_Roe);
      r4[1] = (-0.5*rho_Roe/a_Roe) * ( u_Roe - a_Roe*nxhat(j,i) );
      r4[2] = (-0.5*rho_Roe/a_Roe) * ( v_Roe - a_Roe*nyhat(j,i) );
      r4[3] = (-0.5*rho_Roe/a_Roe) * ( ht_Roe - a_Roe*U_hat_Roe );

      // compute delta's
      drho = rho_R(j,i) - rho_L(j,i);
      du = u_R(j,i) - u_L(j,i);
      dv = v_R(j,i) - v_L(j,i);
      dp = p_R(j,i) - p_L(j,i);

      // compute characteristic variables
      dw1 = drho - dp/(a_Roe*a_Roe);
      dw2 = du*nyhat(j,i) - dv*nxhat(j,i);
      dw3 = du*nxhat(j,i) + dv*nyhat(j,i) + dp/(rho_Roe*a_Roe);
      dw4 = du*nxhat(j,i) + dv*nyhat(j,i) - dp/(rho_Roe*a_Roe);

      U_hat_Roe_L = ( u_L(j,i)*nxhat(j,i) + v_L(j,i)*nyhat(j,i));
      U_hat_Roe_R = ( u_R(j,i)*nxhat(j,i) + v_R(j,i)*nyhat(j,i));

      // grab the first order contributions
      F_L[0] = rho_L(j,i)*U_hat_Roe_L;
      F_R[0] = rho_R(j,i)*U_hat_Roe_R;

      F_L[1] = rho_L(j,i)*u_L(j,i)*U_hat_Roe_L + p_L(j,i)*nxhat(j,i);
      F_R[1] = rho_R(j,i)*u_R(j,i)*U_hat_Roe_R + p_R(j,i)*nxhat(j,i);

      F_L[2] = rho_L(j,i)*v_L(j,i)*U_hat_Roe_L + p_L(j,i)*nxhat(j,i);
      F_R[2] = rho_R(j,i)*v_R(j,i)*U_hat_Roe_R + p_R(j,i)*nxhat(j,i);

      F_L[3] = rho_L(j,i)*ht_L*U_hat_Roe_L;
      F_R[3] = rho_R(j,i)*ht_R*U_hat_Roe_R;

      // Combine 1st order and 2nd order


      // sum over the right vectors

      for (int eq = 0; eq < NEQ; eq++)
      {
        sum2ndOrder=0.5*(abs(lambda1)*dw1*r1[eq]
            + abs(lambda2)*dw2*r2[eq]
            + abs(lambda3)*dw3*r3[eq]
            + abs(lambda4)*dw4*r4[eq]);
        Flux[eq](j,i) = 0.5*(F_L[eq]+F_R[eq]) - sum2ndOrder;
      }
    }
  }
}

void compute2DFlux(
    // This function will take in the Left Right Bottom and Top States,
    // and compute F flux, and G flux respectively.
    MatrixXd* F,           // output - F dir flux (x)         
    MatrixXd* G,           // output - G dir flux (y)
    MatrixXd& n_i_xhat,    // input - xhat normal in i
    MatrixXd& n_i_yhat,    // input - yhat normal in i
    MatrixXd& n_j_xhat,    // input - xhat normal in j
    MatrixXd& n_j_yhat,    // input - yhat normal in j
    MatrixXd* V_L,         // input - Left state prim var
    MatrixXd* V_R,         // input - Right state prim var
    MatrixXd* V_B,         // input - Bottom state prim var
    MatrixXd* V_T,         // input - Top state prim var
    constants C            // input - constants for Roe or VL flags
    )
{

  if (C.f_upwind == 1)
  {
    // i dir
    computeFluxVL(F, n_i_xhat, n_i_yhat, V_L, V_R, C);
    // j dir
    computeFluxVL(G, n_j_xhat, n_j_yhat, V_B, V_T, C);
  }
  else if (C.f_upwind == 2)
  {
    computeFluxRoe(F, n_i_xhat, n_i_yhat, V_L, V_R, C);
    computeFluxRoe(G, n_j_xhat, n_j_yhat, V_B, V_T, C);
  }
  else
  {
    cerr << "ERROR: C.f_upwind Not Defined Correctly!!" << endl;
    exit(1);
  }
}

void applyLimiter(
    // This function avoids repeated code chucks it takes in the iteration number and output 
    // Psi's
    double& psi_p,         // output - modify psi_p double value
    double& psi_n,         // output - modify psi_n double value
    double r_p,            // input - r_p slope 
    double r_n,            // input - r_n slope
    constants C            // input - constants for limiter flag
    ) 
{
  switch (C.f_limiter)
  {
    case 0:  // no limiter
      psi_p = 1.0;
      psi_n = 1.0;
      break;
    case 1: // van leer
      psi_p = (r_p + abs(r_p)) / (1 + abs(r_p));
      psi_n = (r_n + abs(r_n)) / (1 + abs(r_n));
      break;
    case 2: // van albaada
      psi_p = (r_p + r_p*r_p) / (1 + r_p*r_p);
      psi_n = (r_n + r_n*r_n) / (1 + r_n*r_n);
      break;
    case 3: // Ospre
      psi_p = 1.5*(r_p*r_p + r_p) / (1 + r_p + r_p*r_p);
      psi_n = 1.5*(r_n*r_n + r_n) / (1 + r_n + r_n*r_n);
      break;
    case 4: // monotized central least diffusive
      psi_p = mymax(0,mymin(2.0*r_p,mymin(0.5*(1+r_p),2)));
      psi_n = mymax(0,mymin(2.0*r_n,mymin(0.5*(1+r_n),2)));
      break;
    case 5: // MINMOD
      psi_p = mymax(0,mymin(1.0,r_p));
      psi_n = mymax(0,mymin(1.0,r_n));
      break;
    case 6: // superbee
      psi_p = mymax(0,mymax(mymin(2.0*r_p,1.0),mymin(r_p,2.0)));
      psi_n = mymax(0,mymax(mymin(2.0*r_n,1.0),mymin(r_n,2.0)));
      break;
  }
}

void computeUpwindVBT(
    // This function will compute the top and bot state through 
    // similary method as VLR 
    MatrixXd* V_B,         // output - Bottom state prim var   
    MatrixXd* V_T,         // output - Top state prim var
    MatrixXd* V,           // input - prim var
    constants C            // input - constants
    )
{
  int ni_g = V[rhoid].cols();
  int nj_g = V[rhoid].rows(); 
  int ni = ni_g - 2*C.num_ghost;
  int nj = nj_g - 2*C.num_ghost;

  // create the psi pos and neg
  MatrixXd* Psi_Pos = new MatrixXd[NEQ];
  MatrixXd* Psi_Neg = new MatrixXd[NEQ];
  for (int eq = 0; eq < NEQ; eq++)
  {
    Psi_Pos[eq].resize(nj_g+1,ni_g);// interfaces in j
    Psi_Neg[eq].resize(nj_g+1,ni_g);
  }

  // Calculate Psi's  then update top and bottom

  double delta = 1.0e-6;
  for (int i = 0; i < ni_g; i++)
  {
    for (int j = C.num_ghost-2 ; j <= nj+C.num_ghost; j++)
    { 
      // Loop over all of the primvars
      for (int eq = 0; eq < NEQ; eq++)
      {
        // NOW SWEEPING IN THE i direction
        // Calculate the denom which is the same for pos and neg
        double denom = V[eq](j+1,i) - V[eq](j,i);
        // Ensure not division by 0
        denom = SIGN( mymax(abs(denom), delta), denom);
        // Calculate the slopes (r) for both pos and neg dir
        double r_p = (V[eq](j+2,i) - V[eq](j+1, i))/denom;
        double r_n = (V[eq](j,i) - V[eq](j-1, i))/denom;
        // Apply Limiters
        double psi_p, psi_n;
        applyLimiter(psi_p, psi_n, r_p, r_n, C);
        Psi_Pos[eq](j+1,i) = psi_p;
        Psi_Neg[eq](j+1,i) = psi_n;
      }
    }
  }

  for (int i = 0; i < V_T[rhoid].cols(); i++)
  {
    for (int j = 0; j < V_T[rhoid].rows(); j++)
    {
      int j_int = C.num_ghost + j;
      int j_cells = (C.num_ghost - 1) + j;
      int i_cells = C.num_ghost + i;
      for (int eq = 0 ; eq < NEQ; eq++)
      {
        // Note i_cells + 1 is the index for the i+1/2 value of psi
        // Left state at the interface
        V_B[eq](j,i) = V[eq](j_cells,i_cells) + 
          0.25*C.f_eps*((1.0-C.f_kap)*Psi_Pos[eq](j_int-1,i_cells)*
              (V[eq](j_cells,i_cells)-V[eq](j_cells-1,i_cells)) + 
              (1.0+C.f_kap)*Psi_Neg[eq](j_int,i_cells)*
              (V[eq](j_cells+1,i_cells)-V[eq](j_cells,i_cells)));
        // Right state at the interface
        V_T[eq](j,i) = V[eq](j_cells+1,i_cells) - 
          0.25*C.f_eps*((1.0-C.f_kap)*Psi_Neg[eq](j_int+1,i_cells)*
              (V[eq](j_cells+2,i_cells)-V[eq](j_cells+1,i_cells)) 
              + (1.0+C.f_kap)*Psi_Pos[eq](j_int,i_cells)*
              (V[eq](j_cells+1,i_cells)-V[eq](j_cells,i_cells)));
      }  
    }
  }
}

void computeUpwindVLR(
    // Apply upwinding scheme in the i direction to find the left and
    // right states, note it does not use freezing
    MatrixXd* V_L,         // output - Left state prim var     
    MatrixXd* V_R,         // output - right state prim var
    MatrixXd* V,           // input - prim var
    constants C            // input - constants for limiter
    )
{
  int ni_g = V[rhoid].cols();
  int nj_g = V[rhoid].rows(); 
  int ni = ni_g - 2*C.num_ghost;
  int nj = nj_g - 2*C.num_ghost;

  // create the psi pos and neg
  MatrixXd* Psi_Pos = new MatrixXd[NEQ];
  MatrixXd* Psi_Neg = new MatrixXd[NEQ];
  for (int eq = 0; eq < NEQ; eq++)
  {
    Psi_Pos[eq].resize(nj_g,ni_g+1);
    Psi_Neg[eq].resize(nj_g,ni_g+1);
  }


  // Calculate the Psi_Pos and Psi_Neg then update
  double delta = 1.0e-6; // avoid div by zero
  for (int j = 0; j < nj_g; j++)
  {
    for (int i = C.num_ghost-2 ; i <= ni+C.num_ghost; i++)
    {
      //cout << i << endl;
      // Loop over all of the primvars
      for (int eq = 0; eq < NEQ; eq++)
      {
        // Calculate the denom which is the same for pos and neg
        double denom = V[eq](j,i+1) - V[eq](j,i);
        // Ensure not division by 0
        denom = SIGN( mymax(abs(denom), delta), denom);
        // Calculate the slopes (r) for both pos and neg dir
        double r_p = (V[eq](j,i+2) - V[eq](j, i+1))/denom;
        double r_n = (V[eq](j,i) - V[eq](j, i-1))/denom;
        // Apply Limiters
        double psi_p, psi_n;
        applyLimiter(psi_p, psi_n, r_p, r_n, C);
        Psi_Pos[eq](j,i+1) = psi_p;
        Psi_Neg[eq](j,i+1) = psi_n;
      }
    }
  }

  // Fill in the Left and right state values
  for (int i = 0; i < V_L[rhoid].cols(); i++)
  {
    for (int j = 0; j < V_L[rhoid].rows(); j++)
    {
      int i_int = C.num_ghost + i; // loop over i interfaces 
      int i_cells = (C.num_ghost - 1) + i; // loop over j interfaces
      int j_cells = C.num_ghost+j;
      for (int eq = 0 ; eq < NEQ; eq++)
      {
        // Note i_cells + 1 is the index for the i+1/2 value of psi
        // Left state at the interface
        V_L[eq](j,i) = V[eq](j_cells,i_cells) + 
          0.25*C.f_eps*((1.0-C.f_kap)*Psi_Pos[eq](j_cells,i_int-1)*
              (V[eq](j_cells,i_cells)-V[eq](j_cells,i_cells-1)) + 
              (1.0+C.f_kap)*Psi_Neg[eq](j_cells,i_int)*
              (V[eq](j_cells,i_cells+1)-V[eq](j_cells,i_cells)));
        // Right state at the interface
        V_R[eq](j,i) = V[eq](j_cells,i_cells+1) - 
          0.25*C.f_eps*((1.0-C.f_kap)*Psi_Neg[eq](j_cells,i_int+1)*
              (V[eq](j_cells,i_cells+2)-V[eq](j_cells,i_cells+1)) 
              + (1.0+C.f_kap)*Psi_Pos[eq](j_cells,i_int)*
              (V[eq](j_cells,i_cells+1)-V[eq](j_cells,i_cells)));
      }  
    }
  }
}


void MUSCL(
    // This function will fill in the values of the LRBT states used in 
    // the flux calculation it will require 
    // it does not use limiter freezing
    MatrixXd* V_L,         // output - Left state prim var     
    MatrixXd* V_R,         // output - right state prim var
    MatrixXd* V_B,         // output - bottom state prim var
    MatrixXd* V_T,         // output - top state prim var
    MatrixXd* V,           // input - prim var at cells w/g
    constants C            // input - constants for flags
    )
{
  computeUpwindVLR(V_L, V_R, V, C);
  computeUpwindVBT(V_B, V_T, V, C);
}

void setBCMMS(
    // This function just applying the boundary condition of the MMS case
    // It basically copies the solution to the ghost cells
    MatrixXd* V,           // output - Solution at the Boundaries MMS   
    MatrixXd* V_MMS,       // input - Solution at bndries analytically
    constants C            // input - constants for?
    )
{
  int ni_g = V[rhoid].cols(); 
  int nj_g = V[rhoid].rows(); // cell number w/ ghost
  int ni = ni_g - 2*C.num_ghost;
  int nj = nj_g - 2*C.num_ghost;// cell number w/o ghost

  int bot = C.num_ghost - 1; // 1st ghost cell based off 
  int top = C.num_ghost + nj; // top index start
  int left = C.num_ghost -1; // same as bot
  int right = C.num_ghost + ni; // 1st ghost on right

  // j line push in i dir 
  for (int i = C.num_ghost; i < ni+C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      for (int eq = 0; eq < NEQ; eq++)
      {
        // bottom index - j is the filling of ghost
        V[eq](bot-j,i) = V_MMS[eq](bot-j,i);
        // top index top + j is filling ghost
        V[eq](top+j,i) = V_MMS[eq](top+j, i);
      }
    }
  }
  // i line push in j dir 
  for (int i = 0; i < C.num_ghost; i ++)
  {
    for (int j = C.num_ghost; j < C.num_ghost+nj; j++)
    {
      for (int eq = 0; eq < NEQ; eq++)
      {

        // left index - i is the filling of ghost
        V[eq](j,left-i) = V_MMS[eq](j,left-i);
        // right index + i is the filling of ghost
        V[eq](j,right+i) = V_MMS[eq](j,right+i);
      }
    }
  }
}

void solveSolutionMMS(
    // This function evaluates the KNOWN solution which can be used 
    // in determining error
    MatrixXd* V_MMS,       // output - MMS solution at all cells  
    MatrixXd& xc_g,        // input - xcoord cells w/ ghost
    MatrixXd& yc_g,        // input - ycoord cells w/ghost
    constants C            // input - constants for?
    )
{
  // Setup vars associated with Appendix A from the Paper
  double rho0, rhox, rhoy;
  double uvel0, uvelx, uvely;
  double vvel0, vvelx, vvely;
  double wvel0, wvelx, wvely;
  double press0, pressx, pressy;
  double Pi = PI;// Used for mathematica notebook notation
  double gamma = GAMMA;// used for nb notation

  double x, y; // x and y for the sin and cos funcs
  double L = 1.0; // length specified from paper for m.nb

  int ni = xc_g.cols(); // number of cells i dir
  int nj = xc_g.rows(); // number of cells j dir

  // Check for dimension mismatch
  if (V_MMS[rhoid].rows() != nj || V_MMS[rhoid].cols() != ni)
  {
    cerr << "ERROR: Dimension Mismatch in MMS Source Setup?!?!" << endl;
    exit(1);
  }
  if (C.f_supersonic == 0) // subsonic table A.2
  {
    rho0=1.0; rhox=0.15; rhoy=-0.1; // same in both
    uvel0=70.0; uvelx=5.0; uvely=-7.0;
    vvel0=90.0; vvelx=-15.0; vvely=8.5;
    wvel0=0.0; wvelx=0.0; wvely=0.0;
    press0=1.0e5; pressx=0.2e5; pressy=0.5e5; // same in both
  }
  else if (C.f_supersonic == 1) // supersonic table A.1
  {
    rho0=1.0; rhox=0.15; rhoy=-0.1;
    uvel0=800.0; uvelx=50.0; uvely=-30.0;
    vvel0=800.0; vvelx=-75.0; vvely=40.0;
    wvel0=0.0; wvelx=0.0; wvely=0.0;
    press0=1.0e5; pressx=0.2e5; pressy=0.5e5;
  }
  else
  {
    cerr << "ERROR: C.f_supersonic Incorrect Value MMS Source Setup?!?!" << endl;
    exit(1);
  }
  // loop over coords
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      // grab the coords for copying purposes from mathematica
      x = xc_g(j,i); y = yc_g(j,i);

      V_MMS[rhoid](j,i) = rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L);
      V_MMS[uid](j,i) = uvel0 + uvely*cos((3.*Pi*y)/(5.*L)) + uvelx*sin((3.*Pi*x)/(2.*L));

      V_MMS[vid](j,i) = vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2.*Pi*y)/(3.*L));
      V_MMS[pid](j,i) = press0 + pressx*cos((2.*Pi*x)/L) + pressy*sin((Pi*y)/L);
    }
  }
}

void solveSourceMMS(
    // This function will determine the source term that is found from 
    // substituting the MMS solution into the PDE. 
    MatrixXd* S,           // output - Source Term Values                  
    MatrixXd& xc,          // input - x cell coords w/o ghost
    MatrixXd& yc,          // input - y cell coords w/o ghost
    constants C            // input - constants for flags
    )
{
  // Setup vars associated with Appendix A from the Paper
  double rho0, rhox, rhoy;
  double uvel0, uvelx, uvely;
  double vvel0, vvelx, vvely;
  double wvel0, wvelx, wvely;
  double press0, pressx, pressy;
  double Pi = PI;// Used for mathematica notebook notation
  double gamma = GAMMA;// used for nb notation

  double x, y; // x and y for the sin and cos funcs
  double L = 1.0; // length specified from paper for m.nb

  int ni = xc.cols(); // number of cells i dir
  int nj = xc.rows(); // number of cells j dir

  // Check for dimension mismatch
  if (S[rhoid].rows() != nj || S[rhoid].cols() != ni)
  {
    cerr << "ERROR: Dimension Mismatch in MMS Source Setup?!?!" << endl;
    exit(1);
  }
  if (C.f_supersonic == 0) // subsonic table A.2
  {
    rho0=1.0; rhox=0.15; rhoy=-0.1; // same in both
    uvel0=70.0; uvelx=5.0; uvely=-7.0;
    vvel0=90.0; vvelx=-15.0; vvely=8.5;
    wvel0=0.0; wvelx=0.0; wvely=0.0;
    press0=1.0e5; pressx=0.2e5; pressy=0.5e5; // same in both
  }
  else if (C.f_supersonic == 1) // supersonic table A.1
  {
    rho0=1.0; rhox=0.15; rhoy=-0.1;
    uvel0=800.0; uvelx=50.0; uvely=-30.0;
    vvel0=800.0; vvelx=-75.0; vvely=40.0;
    wvel0=0.0; wvelx=0.0; wvely=0.0;
    press0=1.0e5; pressx=0.2e5; pressy=0.5e5;
  }
  else
  {
    cerr << "ERROR: C.f_supersonic Incorrect Value MMS Source Setup?!?!" << endl;
    exit(1);
  }

  // APPLY MMS from the mathematica notebook. Note the same usage of variable names
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      // In the notebook it uses x and y so grab vals
      x = xc(j,i); y = yc(j,i);

      S[rhoid](j,i) = (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))/(2.*L) + (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))/(3.*L) + (Pi*rhox*cos((Pi*x)/L)*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L - (Pi*rhoy*sin((Pi*y)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*
                L)) + vvely*sin((2*Pi*y)/(3.*L))))/(2.*L);

      S[rhouid](j,i) = (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L + (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/(3.*L) + (Pi*rhox*cos((Pi*x)/L)*pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2))/L - (2*Pi*pressx*sin((2*Pi*x)/L))/L - (Pi*rhoy*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*sin((Pi*y)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/ (2.*L) - (3*Pi*uvely*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*sin((3*Pi*y)/(5.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(5.*L);

      S[rhovid](j,i) = (Pi*pressy*cos((Pi*y)/L))/L - (Pi*vvelx*sin((Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/(2.*L) + (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(2.*L) + (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(3.*L) + (Pi*rhox*cos((Pi*x)/L)*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/L - (Pi*rhoy*sin((Pi*y)/(2.*L))*pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/(2.*L); 

      S[rhoetid](j,i) = (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*((-2*Pi*pressx*sin((2*Pi*x)/L))/L + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*((-2*Pi*pressx*sin((2*Pi*x)/L))/((-1 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))) + ((3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L))))/L - (Pi*vvelx*sin((Pi*x)/(2.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/L)/2. - (Pi*rhox*cos((Pi*x)/L)*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L)))/((-1 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L),2))) + (Pi*rhox*cos((Pi*x)/L)*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))))/L) + (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L) + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))))))/(2.*L) + (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L) + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))))))/(3.*L) + (vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)))*((Pi*pressy*cos((Pi*y)/L))/L - (Pi*rhoy*sin((Pi*y)/(2.*L))*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L))/((-1 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L)))))/(2.*L) + (rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))*((Pi*pressy*cos((Pi*y)/L))/((-1 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L))) + ((-6*Pi*uvely*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) + uvelx*sin((3*Pi*x)/(2.*L)))*sin((3*Pi*y)/(5.*L)))/(5.*L) + (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*(vvel0 + vvelx*cos((Pi*x)/(2.*L)) + vvely*sin((2*Pi*y)/(3.*L))))/(3.*L))/2. + (Pi*rhoy*sin((Pi*y)/(2.*L))*(press0 + pressx*cos((2*Pi*x)/L) + pressy*sin((Pi*y)/L)))/(2.*(-1 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.*L)) + rhox*sin((Pi*x)/L),2)))); 
    }
  }
  // FROM the paper the initial values
}

void outputArray(
    // This function will output any eigen matrix into a file
    string FileName,       // input - File
    MatrixXd& out,         //input - matrix
    int n                  //input - iteration
    )
{
  ostringstream StrConvert;
  StrConvert << n; // import n as an ostringstream
  string num_iter = StrConvert.str();// convert to string
  string Address = "./debug/";// add together
  string Suffix = ".txt";
  string Sub = "_";
  Address = Address + FileName + Sub + num_iter + Suffix; // can combine string
  ofstream outfile; // output
  outfile.open(Address.c_str()); // access string component

  // Alot of options for outputting this is just nice to see 
  // where lines end and start
  IOFormat CleanFmt(14,0,", ", "\n", "[", "]");
  outfile << out.transpose().format(CleanFmt) << endl; // output trans
  //outfile << setprecision(14) << out << endl; // output trans
}


void computeTemperature(
    // This function grabs the pressure and dens and computes temperature
    // It is helpful when applying boundary conditions
    MatrixXd& T,           // output - Temperature in K for cells w/ g
    MatrixXd* V            // input - prim var cells w/ g
    )
{
  for (int i = 0; i < T.cols(); i++)
    for (int j =0; j < T.rows(); j++)
      T(j,i) = V[pid](j,i)/(R*V[rhoid](j,i));
}

void initialize(
    // This function initializes all of the primitive variable and it c
    // hecks for the correct case to initialize with. 
    MatrixXd* V,           // output - Primitive variables 
    MatrixXd& xc_g,        // input - xcoord with ghosts
    MatrixXd& yc_g,        // input - ycoord with ghosts
    constants C            // input - constants for case number
    )
{
  int ni = xc_g.cols(); // n cells in i dir
  int nj = xc_g.rows();  // n cells in j dir
  // check for the same sizes
  if (V[rhoid].cols() != ni || V[rhoid].rows() != nj)
  {
    cerr << "ERROR: Dimension Mismatch in Initialization?!?!" << endl;
    exit(1);
  }
  // create important constants specified from papers
  double M;
  double rho;
  double p;
  double u;
  double v;
  double T;
  double AOA;
  // Use a c++ switch for the three cases to save a bunch of if statements
  switch (C.f_case)
  {
    case 1: /////////////// Curvilinear /////////////////
      if (C.f_supersonic == 0) // subsonic
      {
        // Taken from Table A.1 & A.2 from appendix B
        rho = 1.0; // kg/m^3
        u = 70.0; // m/s
        v = 90.0; // m/s
        p = 1.0e5; // Pa
      }
      else if (C.f_supersonic == 1) // supersonic
      {
        rho = 1.0; // kg/m^3
        u = 800.0;  // m/s
        v = 800.0; // m/s
        p = 1.0e5; // Pa
      }
      break;
    case 2: ////////////// 30 Degree Inlet //////////////
      M = 4.0; // na
      T = 217; // K
      // set vars from project desc sheet
      p = 12270.0; // Pa
      rho = p/(R*T); // kg/m^3
      u = M*sqrt(GAMMA*R*T); // m/s
      v = 0.0;
      break;
    case 3: /////////////// NACA Airfoil //////////////
      if (C.f_AOA == 1) // 0 deg
      {
        M = 0.84;
        p = 65855.8; // Pa
        AOA = 0; // deg
      }
      else if (C.f_AOA == 2) // 8 deg
      {
        M = 0.75;
        p = 67243.5; // Pa
        AOA = 8; // deg
      }
      // set vars
      T = 300; // K
      rho = p/(R*T);
      u = M*sqrt(GAMMA*R*T)*cos(AOA*PI/180.0);
      v = M*sqrt(GAMMA*R*T)*sin(AOA*PI/180.0);
      break;
  }
  // fill in the prim var
  V[rhoid] = MatrixXd::Constant(nj,ni,rho); // sets all values to rho
  V[uid] = MatrixXd::Constant(nj,ni,u); // sets all values to u
  V[vid] = MatrixXd::Constant(nj,ni,v); // sets all values to v
  V[pid] = MatrixXd::Constant(nj,ni,p); // sets all values to p
}

void extrapCopyCoords(
    // This function will copy the interior coord of the cells, but
    // will also extrapolate those coords to the ghost cells
    MatrixXd& xc_g,        // output - cell x coord with ghosts 
    MatrixXd& yc_g,        // output - cell y coord with ghosts
    MatrixXd& xc,          // input - cell x coord without ghosts
    MatrixXd& yc,          // input - cell y coord without ghosts
    constants C            // input - constants C for num ghost
    )
{
  int ni = xc.cols();
  int nj = xc.rows();
  // copy the interior values for both x and y
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      int i_g = i+C.num_ghost; // set iter to 3+i
      int j_g = j+C.num_ghost; // set iter to 3+j
      xc_g(j_g, i_g) = xc(j,i);
      yc_g(j_g, i_g) = yc(j,i);
    }
  }
  int bot = C.num_ghost - 1; // 1st ghost cell based off 
  int top = C.num_ghost + nj; // top index start
  int left = C.num_ghost -1; // same as bot
  int right = C.num_ghost + ni; // 1st ghost on right
  // extrapolate to the ghost cells now

  // j dir which means the interior need extrap
  for (int i = C.num_ghost; i < ni+C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bottom index - j is the filling of ghost
      xc_g(bot-j,i) = 2.0*xc_g(bot-j+1,i) - xc_g(bot-j+2,i);
      yc_g(bot-j,i) = 2.0*yc_g(bot-j+1,i) - yc_g(bot-j+2,i);
      // top index top + j is filling ghost
      xc_g(top+j,i) = 2.0*xc_g(top+j-1,i) - xc_g(top+j-2,i);
      yc_g(top+j,i) = 2.0*yc_g(top+j-1,i) - yc_g(top+j-2,i);
    }
  }
  // i dir which means the interior need extrap
  for (int i = 0; i < C.num_ghost; i ++)
  {
    for (int j = C.num_ghost; j < C.num_ghost+nj; j++)
    {
      // left index - i is the filling of ghost
      xc_g(j,left-i) = 2.0*xc_g(j,left-i+1) - xc_g(j,left-i+2);
      yc_g(j,left-i) = 2.0*yc_g(j,left-i+1) - yc_g(j,left-i+2);
      // right index + i is the filling of ghost
      xc_g(j,right+i) = 2.0*xc_g(j,right+i-1) - xc_g(j,right+i-2);
      yc_g(j,right+i) = 2.0*yc_g(j,right+i-1) - yc_g(j,right+i-2);
    }
  }
  //cout << xc_g << endl;
  //set to 1 out the corners for sake of clarity.
  for (int i = 0; i < C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bot left
      xc_g(j,i) = 1.0; yc_g(j,i) = 1.0;
      // bot right
      xc_g(j,right+i) = 1.0; yc_g(j,right+i) = 1.0;
      // top right
      xc_g(top+j,right+i) = 1.0; yc_g(top+j,right+i) = 1.0;
      // top left
      xc_g(top+j,i) = 1.0; yc_g(top+j,i) = 1.0;
    }
  }
  //cout << xc_g << endl;
}
void computeVolume(
    // This function computes the volumes of the interior cell points. 
    MatrixXd& Volume,      // output - Volume of the cells
    MatrixXd& xn,          // input - x coord for nodal vals
    MatrixXd& yn           // input - y coord for nodal vals
    )
{
  double dx1, dx2, dy1, dy2; // find the diagonals on the cell
  double cross; // cross product
  for (int i = 0; i < Volume.cols(); i++)
  {
    for (int j = 0; j < Volume.rows(); j++)
    {
      dx1 = xn(j,i) - xn(j+1,i+1); // bot L - top R
      dy1 = yn(j,i) - yn(j+1,i+1); // bot L - top R
      dx2 = xn(j+1,i) - xn(j,i+1); // top L - bot R
      dy2 = yn(j+1,i) - yn(j,i+1); // top L - bot R
      cross = abs(dx1*dy2 - dx2*dy1); // cross product
      Volume(j,i) = 0.5*cross; // assume w=1 volume of cell
    }
  } 
}

void computeNormalVectors(
    // This function takes in the previously calculated areas and
    // Finds the fraction in the x and y physical directions
    MatrixXd& n_i_xhat,    // output - i with A of normal vector in x phys
    MatrixXd& n_i_yhat,    // output - i comp but the yhat phys
    MatrixXd& n_j_xhat,    // output - j with A of normal vector in x phys
    MatrixXd& n_j_yhat,    // output - j dir but with yhat phys
    MatrixXd& xn,          // input - nodal x coordinate
    MatrixXd& yn,          // input - nodal y coordinate
    MatrixXd& Ai,          // input - Area pointed in i (x square grid)
    MatrixXd& Aj           // input - Area pointed in j (y square grid)
    )
{
  // NOTE i is the columns and j are the rows for A(row,col)
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      // See notes Sec 6 slide 75
      n_i_xhat(j,i) = (yn(j+1,i) - yn(j,i)) / Ai(j,i); 
      n_i_yhat(j,i) = -(xn(j+1,i) - xn(j,i)) / Ai(j,i);
    }
  }
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      n_j_xhat(j,i) = -(yn(j,i+1) - yn(j,i)) / Aj(j,i);
      n_j_yhat(j,i) = (xn(j,i+1) - xn(j,i)) / Aj(j,i);
    }
  }
  //cout << n_i_yhat << endl;
}

void computeArea(
    // This function computes the area based on the node locations 
    // for both the i facing direction and the j facing direction
    MatrixXd& Ai,          // output - i-direction areas (c, n)size    
    MatrixXd& Aj,          // output - j-direction areas (n, c)size
    MatrixXd& xn,          // input - x coord for nodes 
    MatrixXd& yn           // input - y coord for nodes
    )
{
  // loop for Ai
  // NOTE i is the columns and j is the rows
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      double dx = xn(j,i) - xn(j+1,i); // compute bottom - top x
      double dy = yn(j,i) - yn(j+1,i); // compute bottom - top y
      Ai(j,i) = sqrt( dx*dx + dy*dy ); // pythag
    }
  }
  // loop for Aj
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      double dx = xn(j,i) - xn(j,i+1); // compute left - right x
      double dy = yn(j,i) - yn(j,i+1); // compute left - right y
      Aj(j,i) = sqrt( dx*dx + dy*dy);  // pytha
    }
  }
  // NOTE area is the dist times a w into the z dir which is just 1 for 2d
}




void inputMesh(
    // This function will fill in the x, y, and z values read from the grid files and return
    // nodal coords and also cell center coords 
    MatrixXd& xn,          // output - xcoord for all nodes 
    MatrixXd& yn,          // output - ycoord for all nodes
    MatrixXd& zn,          // output - zcoord for all nodes
    MatrixXd& xc,          // output - xcoord for all centers
    MatrixXd& yc,          // output - ycoord for all centers
    MatrixXd& zc,          // output - zcoord for all centers
    constants C            // input - Constants for determining which grid file
    )
{
  double readval; // read in value
  vector<double> data; // create an abitrarily long vector
  string line; // string of the line 

  // get filename
  string filename = readMeshName(C);
  const char *FILENAME = filename.c_str();
  // open the instream
  ifstream infile;
  infile.open(FILENAME);

  // check for reading in error
  if (infile.fail()) 
  {
    cerr << "ERROR: Mesh File Cannot Be Opened!\nDid you extract the grids?!?!" << endl;
    exit(1);
  }

  // Read in data from file
  while (infile.good()) 
  {
    while (getline(infile, line)) 
    {
      istringstream streamA(line);
      while (streamA >> readval)
      {
        data.push_back(readval);
      }
    }
  }
  // data[0] is some number thats irrelevant
  int ni = data[1]-1; // 2nd number in the file is the number of cols (xvals i)
  int nj = data[2]-1; // 3rd number in the file is the number of rows (yvals j)
  int nk = data[3]-1; // number of z vals

  // resize output based off of the read input
  xc.resize(nj,ni); // resize (rows, cols) 4 rows by 5 cols : nj=4, ni=5
  yc.resize(nj,ni);
  zc.resize(nj,ni);

  xn.resize(nj+1,ni+1);
  yn.resize(nj+1,ni+1);
  zn.resize(nj+1,ni+1);

  int x_index = 4; // start of x vals
  int y_index = x_index + (nj+1)*(ni+1)*(nk+1);
  int z_index = y_index + (nj+1)*(ni+1)*(nk+1);


  // fill in the nodal values

  for (int row = 0; row < xn.rows(); row++)
  {
    for (int col = 0 ; col < xn.cols(); col++)
    {
      xn(row,col) = data[x_index];
      yn(row,col) = data[y_index];
      zn(row,col) = data[z_index];
      x_index++;
      y_index++;
      z_index++;
    }
  }

  for (int row = 0; row < xc.rows(); row++)
  {
    for (int col = 0 ; col < xc.cols(); col++)
    {
      // take average of all four xvals from corn
      // order: botL + topL + botR + topR 
      xc(row,col) = 0.25*(xn(row,col) + xn(row+1,col) + 
          xn(row,col+1) + xn(row+1,col+1));
      yc(row,col) = 0.25*(yn(row,col) + yn(row+1,col) + 
          yn(row,col+1) + yn(row+1,col+1));
      zc(row,col) = 0.25*(zn(row,col) + zn(row+1,col) + 
          zn(row,col+1) + zn(row+1,col+1));
    }
  }
}

string readMeshName(
    constants C
    )
{

  if (C.f_case == 1) { // curvelinear meshes 
    if (C.f_mesh == 1) // read in smallest mesh
      return "grids/Curvilinear/2d9.grd";
    if (C.f_mesh == 2) 
      return "grids/Curvilinear/2d17.grd";
    if (C.f_mesh == 3) 
      return "grids/Curvilinear/2d33.grd";
    if (C.f_mesh == 4)
      return "grids/Curvilinear/2d65.grd";
    if (C.f_mesh == 5) 
      return "grids/Curvilinear/2d129.grd";
    if (C.f_mesh == 6)
      return "grids/Curvilinear/2d257.grd";
  }
  if (C.f_case == 2) { // 30 deg inlet
    if (C.f_mesh == 1) // read in smallest mesh
      return "grids/Inlet/53x17.grd";
    if (C.f_mesh == 2) 
      return "grids/Inlet/105x33.grd";
    if (C.f_mesh == 3) 
      return "grids/Inlet/209x65.grd";
    if (C.f_mesh == 4)
      return "grids/Inlet/417x129.grd";
  }
  if (C.f_case == 3) { // NACA airfoild
    if (C.f_mesh == 1) // read in smallest mesh
      return "grids/NACA64A006/49x14.grd";
    if (C.f_mesh == 2) 
      return "grids/NACA64A006/97x27.grd";
    if (C.f_mesh == 3) 
      return "grids/NACA64A006/193x53.grd";
    if (C.f_mesh == 4)
      return "grids/NACA64A006/385x105.grd";
  }
  return "failure"; // return failure because it can not recognize the file
}


string buildCaseFolder(
    // This function will build the folder in the output directory corresponding to 
    // the values given from the input files. 
    constants C            // input - Constants contain evry thing from inputfile
    )
{
  // Create converters
  ostringstream StrConvert1; 
  ostringstream StrConvert2;
  ostringstream StrConvert3;
  ostringstream StrConvert4;
  ostringstream StrConvert5;
  // grab and convert values to strings to write in folder name
  StrConvert1 << C.f_case;
  string Case = StrConvert1.str();
  StrConvert2 << C.f_mesh;
  string Mesh = StrConvert2.str();
  StrConvert3 << C.f_upwind;
  string Flux = StrConvert3.str();
  StrConvert4 << C.f_eps+1; // add 1 because it uses 0 and 1
  string Order = StrConvert4.str();
  StrConvert5 << C.f_AOA;
  string AOA = StrConvert5.str();

  string Output_Folder; // construct the name of the output folder based on C
  if (C.f_case != 3) // no AOA
    Output_Folder  = "Case_"+Case+"_Mesh_"+Mesh+"_Flux_"+Flux+"_order_"+Order;
  if (C.f_case == 3) // no AOA
    Output_Folder  = "Case_"+Case+"_Mesh_"+Mesh+"_Flux_"+Flux+"_AOA_"+AOA+"_order_"+Order;

  // Remove any content from directory for the new simulation
  string wipe_dir = "exec rm -rf output/"+Output_Folder+"/*"; // remove everything
  string make_dir = "exec mkdir -p output/"+Output_Folder; // create directory

  const char *exec1 = wipe_dir.c_str(); // convert to type that system takes
  const char *exec2 = make_dir.c_str();

  int systemRet1 = system(exec1); // run the bash command with exec
  int systemRet2 = system(exec2);  

  if(systemRet1 == -1 || systemRet2 == -1)
  {
    cerr << "ERROR: Could Not Execute mkdir Command!" << endl;
    exit(1);
  }

  exec1 = new char; delete exec1; exec1 = NULL;
  exec2 = new char; delete exec2; exec2 = NULL;
  return Output_Folder; // return the OutputFolder for use in writing out
}

double searchInputFile(
    // This function will search through the FILENAME to find the Var name 
    // and grab the value associated with Var if it is enclose with : and 
    // ;
    string FileName,       // input - filename to search through    
    string Var             // input - variable name from file
    )
{
  const char *FILENAME = FileName.c_str(); // convert string to character C arr
  ifstream in(FILENAME); // create inflow stream from file
  stringstream buffer; // read in everything
  buffer << in.rdbuf(); // read buffer associated with inflow
  string text = buffer.str();

  // define cursor
  size_t pos1 = 0; // represents and unsigned int type returned by sizeof operator
  size_t pos2; // right cursor

  // Store Strings in array
  string InputStr;
  string InputName = Var + ":";

  // Find where in the text string from buffer of infile stream is InputName
  pos1 = text.find(InputName, pos1); // grab the position by name
  pos1 = text.find(":", pos1); // Then grab position of : from name
  pos2 = text.find(";", pos1); // search for semicolon pos2 will enclose the variable
  InputStr = text.substr(pos1+1, (pos2-pos1-1)); // cut out the string where number lies
  FILENAME = new char; delete FILENAME; FILENAME = NULL; // FILENAME delete so not used 
  return atof(InputStr.c_str()); // convert the InputStr to type double 
}

constants loadInputFile(
    // This function will load the inputs into the program from an external file so 
    // that cases can be run in parrallel
    string FileName        // input - filename to read in the parameters 
    )
{
  constants C;
  // Grab the values from the search algorithm looking for f_case: MYNUMBER;
  C.f_case = searchInputFile(FileName, "f_case");
  C.f_mesh = searchInputFile(FileName, "f_mesh");
  C.f_supersonic = searchInputFile(FileName, "f_supersonic");
  C.f_AOA = searchInputFile(FileName, "f_AOA");
  C.f_upwind = searchInputFile(FileName, "f_upwind");
  C.f_limiter = searchInputFile(FileName, "f_limiter");
  C.f_eps = searchInputFile(FileName, "f_eps");
  C.rk_order = searchInputFile(FileName, "rk_order");
  C.nmax = searchInputFile(FileName, "nmax");
  C.wint = searchInputFile(FileName, "wint");
  C.pint = searchInputFile(FileName, "pint");
  C.num_ghost = searchInputFile(FileName, "num_ghost");
  C.localdt = searchInputFile(FileName, "localdt");
  C.tol = searchInputFile(FileName, "tol");
  C.cfl = searchInputFile(FileName, "cfl");
  C.f_kap = searchInputFile(FileName, "f_kap");
  return C;
}

double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}
