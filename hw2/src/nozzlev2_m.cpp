/* 
 * Computational Fluid Dynamics
 * This code will compute the exact isentropic solution to the 1-D axisymmetric converging-diverging supersonic nozzle.
 * Written By: Robert Masti
 * 01/27/2017
 * This is the main file that will that will use a header file containing to function prototypes, structure declarations, symbolic constants, class declaration, etc.. And will use a seperate cpp file for the function declarations. 
 */
#include "nozzlev2.hpp"
FILE *fp1; //output the exact solution
FILE *fp2; //output the supersonic outflow final 
FILE *fp3; //output the subsonic outflow
FILE *fp4; //Residuals norms supersonic
FILE *fp5; //Residuals norms subsonic

int main()
{

  output_file_headers();//This prints out the output file headers which can be modified at the bottom of this file.

  ////////////////// Set local constants for the simulation //////////////////
  //Different from the global constants defined in the header file
  constants consts; //constants is a struct declared in the header file
  consts.p0 = 300.0*1000.0; //Pa
  consts.T0 = 600.0; //Kelvin
  consts.A_t = A_x(0.0);//should be m^2
  consts.tol = 1.0e-12; //convergence tolerance 
  consts.gamma = 1.4; 
  consts.outflow = true; //True is supersonic outflow
  consts.cfl = 0.1; //cfl is 0.1 but can be changed after a certain number of iterations

  //////////////////////// START of SIMULATIONS ////////////////////////////
  //This main file is kept purposely really short and it just tells these two files to run their respective sims (exact solution, supersonic outflow, and TBD subsonic function
  
  isentropic_exact(consts);// Pass the constants to isentropic exact and isentropic 
  quasi1Dnozzle(consts); // Run with supersonic outflow

  consts.pb = 120.00*1000.0; //Pa
  consts.outflow = false; //subsonic outflow condition
  quasi1Dnozzle(consts); // Run with subsonic outflow
  //////////////////////// END of SIMULATIONS ////////////////////////////
 
  cout << "Done see data directory" << endl;// Identify that all simulations have finished

  return 0;
  
}

void quasi1Dnozzle(constants C)
{

  // Setup the number of columns
  int number_of_cells = N + 2*num_ghost_cells;
  // 
  // Initialize Matrices
  MatrixXd xi(1,number_of_cells+1);//xcoord for all cell interfaces
  MatrixXd xc(1,number_of_cells);//xcoord for all cell centers
  set_geometry(xc, xi);//Fills in xc and xi

  //Applies the A_x(x) function to all elements of xcenter the same can
  //be done for xinterface if need be;
  MatrixXd Ac = xc.unaryExpr(&A_x);
  MatrixXd Ai = xi.unaryExpr(&A_x);
  //Similarly this can be done for the mach number at the start only for
  //center cell values
  MatrixXd Mc = xc.unaryExpr(&M_xinitial);

  // The following creates 4 U values on the entire grid of 1xnumbcells
  // As an example to access U2 on node 2,1 it would be U[1](1,0)
  MatrixXd* V=new MatrixXd[4]; //primitive variables
  MatrixXd* U=new MatrixXd[4]; //conservative variables
  MatrixXd* Uinterface=new MatrixXd[4];//conservative var at interface
  MatrixXd* F=new MatrixXd[4]; //F flux 
  MatrixXd* d=new MatrixXd[4]; //artifical d flux
  MatrixXd* S=new MatrixXd[4]; //source term
  MatrixXd* Res=new MatrixXd[4]; //residual term only on domain of interest
  //set size of the array 
  for (int i = 0 ; i < 4 ; i++)
  {
    U[i].resize(1,number_of_cells);
    V[i].resize(1,number_of_cells);
    Uinterface[i].resize(1,N+1);
    F[i].resize(1,N+1);
    d[i].resize(1,N+1);
    Res[i].resize(1,N);
    S[i].resize(1,N);
  }

  //initialize primitive variable extrapolate and then apply B.C's;
  initialize(V, Mc, C);

  ///////////////// MAIN LOOP ////////////////
  double t = 0.0;
  double L2normold1 = 10.0;
  double L2normold2 = 10.0;
  double L2normold3 = 10.0;
  double L2normold4 = 10.0;
  for(int n= 0 ;  n < nmax; n++)
  {
    // Beginning of Iteration Loop
    //extrapolate the interior cells to the ghost cells
    extrapolate_to_ghost(V);

    //Apply Boundary Conditions
    set_boundary_conditions(V, C);
   
    //Find lambdamax at cell centers
    MatrixXd Lambda_mcenter(1, number_of_cells);
    compute_lambda(Lambda_mcenter, V, C);
   
    //After extrap and boundary convert V to U
    primtocons(U, V, C);
   
    //Reconstruct to find Uinterface and lambdainterface
    MatrixXd Lambda_minterface(1, N+1);//ONLY AT INTERFACES OF INTEREST
    reconstruct(Uinterface, Lambda_minterface, U, Lambda_mcenter);

    //Now with reconstructed values compute flux
    compute_F_flux(F, Uinterface, C);

    //Find the artificial viscosity flux d
    artificial_viscosity(d, V, U, Lambda_minterface);

    //compute the S source term on the domain
    compute_source(S, V, xc);

    //NOTE Ai is defined at every interface including ghost cells!!
    //Compute residuals to be used in norms and for completing a step
    compute_residual(Res, S, F, d, Ai);

    //Execute 1 time iteration 
    double timestep;
    iteration(U, timestep, Res, Ac, Lambda_mcenter, C);
    t+=timestep;
    //Update the primitive variables with the new U
    constoprim(V, U, C);
    // End of Iteration Loop

    /////////////// OUTPUT SOLNS & RESIDUALS ///////////////
    MatrixXd L2norm(1, neq);
    MatrixXd L2normrel(1,neq);
    compute_norms(L2norm, Res);
    L2normrel(0,0) = L2norm(0,rhoid)/L2normold1;
    L2normrel(0,1) = L2norm(0,uid)/L2normold2;
    L2normrel(0,2) = L2norm(0,vid)*0.0;
    L2normrel(0,3) = L2norm(0,pid)/L2normold4;
    if (n == 3)
    {
      L2normold1 = L2norm(0,rhoid);
      L2normold2 = L2norm(0,uid);
      L2normold3 = L2norm(0,vid);
      L2normold4 = L2norm(0,pid);
    }
    if (L2normrel(0,rhoid) < C.tol && L2normrel(0,uid) < C.tol &&  L2normrel(0,pid) < C.tol)
    {
      //OUTPUT THE SOLUTION
      cout << "CONVERGENCE CRITERIA MET" << endl;
      //write out the final solution for subsonic and supersonic outflow
      if (C.outflow == true)
      {
        //supersonic outflow
        write_solution(fp2, xc, Ac, V, U, C);
        fprintf(fp4, "%i %e %e %e %e %e\n", n, t, L2normrel(0,rhoid), 
            L2normrel(0,uid), L2normrel(0,vid), L2normrel(0,pid));
        fclose(fp2);// Close the exact solution file
        fclose(fp4);// Close the residual history
      }
      else
      {
        //subsonic outflow
        write_solution(fp3, xc, Ac, V, U, C);
        fprintf(fp5, "%i %e %e %e %e %e\n", n, t, L2normrel(0,rhoid), 
            L2normrel(0,uid), L2normrel(0,vid), L2normrel(0,pid));
        fclose(fp3);// Close the exact solution file
        fclose(fp5);
      }
      break;
    }
    if (n % 1000 == 0)
    {
      //OUTPUT RESIDUAL HISTORY
      if (C.outflow == true)
      {
        //supersonic outflow
        fprintf(fp4, "%i %e %e %e %e %e\n", n, t, L2normrel(0,rhoid), 
            L2normrel(0,uid), L2normrel(0,vid), L2normrel(0,pid));
      }
      else
      {
        //subsonic outflow
        fprintf(fp5, "%i %e %e %e %e %e\n", n, t, L2normrel(0,rhoid), 
            L2normrel(0,uid), L2normrel(0,vid), L2normrel(0,pid));
      }
      cout << "n = " << n <<  " L2normrel: rho= " << L2normrel(0,0) << 
        " rhou= " << L2normrel(0,1) << 
        " rhov= " << L2normrel(0,2) <<
        " rhoet= "<< L2normrel(0,3) << endl;
    }
  }
  //delete Res; Res = NULL;
}


void isentropic_exact(constants C)
{
  int number_of_cells = N+2*num_ghost_cells;
  MatrixXd xc(1, N+2*num_ghost_cells);
  MatrixXd xi(1, N+2*num_ghost_cells+1);
  MatrixXd* U=new MatrixXd[4]; //Conserved variables
  MatrixXd* V=new MatrixXd[4]; //primitive variables
  for (int i = 0 ; i < 4 ; i++)
  {
    U[i].resize(1,number_of_cells);
    V[i].resize(1,number_of_cells);
  }
  set_geometry(xc, xi);//Fills in xc and xi
  MatrixXd Ac = xc.unaryExpr(&A_x);
  MatrixXd Mc = xc.unaryExpr(&M_xinitial);
  exactsol(V, Ac, C); //pass the datastruct and Area val to return all primary variables
  primtocons(U, V, C);

  //OUTPUT EXACT SOLUTION 
  write_solution(fp1, xc, Ac, V, U, C);
  fclose(fp1);// Close the exact solution file
}

void output_file_headers()
{
  fp1 = fopen("data/ss/exactsol.dat", "w");
  fprintf(fp1, "Isentropic Exact Solution: I=%i\n", N);
  fprintf(fp1, "x(m) Area(m^2) rho(kg/m^3) u(m/s) v(m/s) Press(Pa) Mach U1 U2 U3 U4\n");

  fp2 = fopen("data/ss/q1Dnozzle.dat","w");
  fprintf(fp2, "Final Quasi-1D Nozzle Supersonic Solution: I=%i\n",N);
  fprintf(fp2, "x(m) Area(m^2) rho(kg/m^3) u(m/s) v(m/s) Press(kPa) Mach U1 U2 U3 U4\n");

  fp3 = fopen("data/sb/q1Dnozzle.dat","w");
  fprintf(fp3, "Final Quasi-1D Nozzle Subsonic Solution: I=%i\n",N);
  fprintf(fp3, "x(m) Area(m^2) rho(kg/m^3) u(m/s) v(m/s) Press(kPa) Mach U1 U2 U3 U4\n");

  fp4 = fopen("data/ss/q1Dhistory.dat","w");
  fprintf(fp4, "Quasi-1D Iterative Convergence Supersonic\n");
  fprintf(fp4, "Iteration Time(s) L2norm1 L2norm2 L2norm3 L2norm4\n");

  fp5 = fopen("data/sb/q1Dhistory.dat","w");
  fprintf(fp5, "Quasi-1D Iterative Convergence Subsonic\n");
  fprintf(fp5, "Iteration Time(s) L2norm1 L2norm2 L2norm3 L2norm4\n");
}
