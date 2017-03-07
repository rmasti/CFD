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

  //Choose which data you wish to generate:
  //isentropic exact soln
  //isentropic fvm
  //shock fvm

  ////////////////// Set local constants for the simulation //////////////////
  //Different from the global constants defined in the header file
  constants consts; //constants is a struct declared in the header file
  consts.p0 = 300.0*1000.0; //Pa
  consts.T0 = 600.0; //Kelvin
  consts.A_t = A_x(0.0);//should be m^2
  consts.cond = true; //This chooses supersonic outflow conditions on the exact solution
  consts.tol = 1.0e-14; //convergence tolerance 
  consts.gamma = 1.4; 
  consts.outflow = true; //True is supersonic outflow
  consts.cfl = 0.1; //cfl is 0.1 but can be changed after a certain number of iterations

  
  output_file_headers();//This prints out the output file headers which can be modified at the bottom of this file.
  //////////////////////// RUN THE SIMULATIONS ////////////////////////////
  //This main file is kept purposely really short and it just tells these two files to run their respective sims (exact solution, supersonic outflow, and TBD subsonic function
  
  //isentropic_exact(consts);// Pass the constants to isentropic exact and isentropic 
  fclose(fp1);// Close the exact solution file

  ////////// SUPERSONIC ///////////
  quasi1Dnozzle(consts); // Run with supersonic outflow
  ////////// SUBSONIC ///////////
  consts.pb = 120.00*1000.0; //Pa
  consts.outflow = false;
  //quasi1Dnozzle(consts); // Run with subsonic outflow
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
  //Similarly this can be done for the mach number at the start only for
  //center cell values
  MatrixXd Mc = xc.unaryExpr(&M_xinitial);
  //cout << xi << endl;
  //cout << xc << endl;
  //cout << Ac << endl;
  //cout << Mc << endl;

  // The following creates 4 U values on the entire grid of 1xnumbcells
  // As an example to access U2 on node 2,1 it would be U[1](1,0)
  //MatrixXd U[4];
  //MatrixXd V[4];
  MatrixXd* V=new MatrixXd[4];
  MatrixXd* U=new MatrixXd[4];
  MatrixXd* Uinterface=new MatrixXd[4];
  MatrixXd* F=new MatrixXd[4];
  for (int i = 0 ; i < 4 ; i++)
  {
    U[i].resize(1,number_of_cells);
    V[i].resize(1,number_of_cells);
    Uinterface[i].resize(1,N+1);
    F[i].resize(1,N+1);
  }
 
  //initialize primitive variable extrapolate and then apply B.C's;
  initialize(V, Mc, C);
  //cout << V[rhoid].transpose() << endl;
  //cout << V[uid].transpose() << endl;
  //cout << V[vid].transpose() << endl;
  //cout << V[pid].transpose() << endl;
 
  //extrapolate the interior cells to the ghost cells
  extrapolate_to_ghost(V);
  //cout << V[rhoid].transpose() << endl;
  
  //Apply Boundary Conditions
  set_boundary_conditions(V, C);
  //cout << V[rhoid] << endl;
  //cout << V[uid] << endl;
  //cout << V[vid] << endl;
  //cout << V[pid] << endl;

  primtocons(U, V, C);
  //cout << U[rhoid] << endl;
  //cout << U[rhouid] << endl;
  //cout << U[rhovid] << endl;
  //cout << U[rhoetid] << endl;


  //Find lambdamax at cell centers
  MatrixXd Lambda_mcenter(1, number_of_cells);
  compute_lambda(Lambda_mcenter, V, C);
  
  //Reconstruct to find Uinterface and lambdainterface
  MatrixXd Lambda_minterface(1, N+1);//ONLY AT INTERFACES OF INTEREST
  reconstruct(Uinterface, Lambda_minterface, U, Lambda_mcenter);
  //cout << Lambda_mcenter << endl;
  //cout << Lambda_minterface << endl;
  cout << Uinterface[rhoid] << endl;
  cout << Uinterface[rhouid] << endl;
  cout << Uinterface[rhovid] << endl;
  cout << Uinterface[rhoetid] << endl;





  //Now with reconstructed values compute flux
  compute_F_flux(F, Uinterface, C);
  cout << F[frhouid] << endl;
  cout << F[frhouuid] << endl;
  cout << F[frhovid] << endl;
  cout << F[frhouhtid] << endl;



 




  //compute_fluxes(

/*
 vector<double> Res3;
  double time=0.0;

  double dt=0.0;
  ////////////////////////// MAIN TIME LOOP ///////////////////////////
  //choose iterator n to be time evolution and iterate until convergence 
  //or until maximum iterations have been reached
  for (int n=0; n < nmax; n++)
  {
    //initialize a flux vector to be filled only on interior cells
    vector<fluxes> Farr; 
    //extrapolate from domain of interest to the ghost cells
    extrapolate_to_ghost(Uold, C);
    //After extrapolating values set the values at the outermost edge
    set_boundary_cond(Uold, C); //good here and up
    //After getting Uold ready for iteration compute fluxes 
    compute_fluxes(Farr, Uold, C);
    //Main iteration step 

    time += dt;
    iteration_step(Unew, Resarr, dt, Marr, Farr, Uold, X_center, X_interface, C);
    //With Resarr compute norms
    compute_norms(Linfnorm, L1norm, L2norm, Resarr); 
    t.push_back(time);
    nvec.push_back(n);
    Res1.push_back(L2norm[0]/Normold[0]);
    Res2.push_back(L2norm[1]/Normold[1]);
    Res3.push_back(L2norm[2]/Normold[2]);

    //Write out the solution after ever __ number of time steps (currently 5000)
    if (n % 500 == 0){
      cout <<"n = " << n <<  " RES: " << L2norm[0]/Normold[0] <<  " " << L2norm[1]/Normold[1] << " " << L2norm[2]/Normold[2] << endl;
      FILE *ftemp;
      //Create necessary header for tecplot
      std::ostringstream ss;
      if (C.outflow == true) //supersonic 
      {
        ss << "data/ss/q1Dnozzle_" << n << ".dat";
        ftemp = fopen(ss.str().c_str(), "w");
        fprintf(ftemp, "TITLE = \"Supesonic Outflow Exact Solution\"\n");  
      }
      else
      {
        ss << "data/sb/q1Dnozzle_" << n << ".dat";
        ftemp = fopen(ss.str().c_str(), "w");
        fprintf(ftemp, "TITLE = \"Subsonic Outflow Exact Solution\"\n"); 
      }
      fprintf(ftemp, "variables=\"x(m)\"\"Area(m^2)\"\"rho(kg/m^3)\"\"u(m/s)\"\"Press(kPa)\"\"Mach\"\"U1\"\"U2\"\"U3\"\n");
      fprintf(ftemp, "zone T=\"%i\"\n", n);
      fprintf(ftemp , "I=%i\n",N);
      fprintf(ftemp, "DATAPACKING=POINT\n");
      fprintf(ftemp, "DT = (DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )\n");
      //fill in the file
      write_out(ftemp, Aarr, X_center,  Unew, C);
      //close the file
      fclose(ftemp);
    } 
    //fill in the old norms so that these can be tracked for convergence
    if (n==3)
    {
      Normold[0] = L2norm[0];
      Normold[1] = L2norm[1];
      Normold[2] = L2norm[2];
    }
    // convergence criteria
    
    if( L2norm[0]/Normold[0] < C.tol && L2norm[1]/Normold[1] < C.tol && L2norm[2]/Normold[2] < C.tol)
    {
      cout<<"CONVERGENCE CRITERIA MET: in n = "<<n<<" iterations"<<endl;
      if (C.outflow == true) 
      {
        write_out(fp2, Aarr, X_center, Unew, C);
        write_res(fp4, nvec, t, Res1, Res2, Res3);
      }
      else 
      {
        write_out(fp3, Aarr, X_center, Unew, C); 
        write_res(fp5, nvec, t, Res1, Res2, Res3);
      }
      break;
    }
    //update Uold with Unew
    Uold = Unew;
  }
  */
    delete U; U = NULL;
    delete V; V = NULL;
    delete F; F = NULL;
    delete Uinterface; Uinterface = NULL;
//    delete Lambda; Lambda = NULL;
//    delete Lambda_interface; Lambda_interface = NULL;
}
/*
void isentropic_exact(constants C)
{
  vector<double> x(N); //array of x values 
  vector<double> Aarr(N); //array of areas that will be passed to the rootfinder funct.
  vector<prim_exact> prvparr(N); //Another data struct defined in header.
  vector<consvar> Uexact(N);
  for (int i = 0; i<N;i++)
  {
    if (i<N/2) C.cond = true; //use a subsonic initial M guess
    else C.cond = false; //use a supersonic initial M guess
    x[i] = xmin + (dx/2)*(i+1); //Evenly spaced x values
    Aarr[i] = A_x(x[i]); //The area is defined in hw1_f.cpp
    prvparr[i] = exactsol(Aarr[i],C); //pass the datastruct and Area val to return all primary variables
    primvar Vtemp; 
    Vtemp.rho = prvparr[i].rho;
    Vtemp.u = prvparr[i].u;
    Vtemp.p = prvparr[i].p; 
    Uexact[i] = primtocons(Vtemp, C);
  }
  write_out(fp1, Aarr, x, Uexact, C); //good
}
*/
/****************************** FILE HANDLING **********************************/
void output_file_headers()
{
  fp1 = fopen("data/ss/exactsol.dat", "w");
  fprintf(fp1, "TITLE = \"Isentropic Exact Solution\"\n");
  fprintf(fp1, "variables=\"x(m)\"\"Area(m^2)\"\"rho(kg/m^3)\"\"u(m/s)\"\"Press(kPa)\"\"Mach\"\"U1\"\"U2\"\"U3\"\n");
  fprintf(fp1 , "I=%i\n",N);
  fprintf(fp1, "DATAPACKING=POINT\n");
  fprintf(fp1, "DT = (DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )\n");
 
  fp2 = fopen("data/ss/q1Dnozzle_final.dat","w");
  fprintf(fp2, "TITLE = \"Final Quasi-1D Nozzle Supersonic Solution\"\n");
  fprintf(fp2, "variables=\"x(m)\"\"Area(m^2)\"\"rho(kg/m^3)\"\"u(m/s)\"\"Press(kPa)\"\"Mach\"\"U1\"\"U2\"\"U3\"\n");
  fprintf(fp2 , "I=%i\n",N);
  fprintf(fp2, "DATAPACKING=POINT\n");
  fprintf(fp2, "DT = (DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )\n");

  fp3 = fopen("data/sb/q1Dnozzle_ss_final.dat","w");
  fprintf(fp3, "TITLE = \"Final Quasi-1D Nozzle Subsonic Solution\"\n");
  fprintf(fp3, "variables=\"x(m)\"\"Area(m^2)\"\"rho(kg/m^3)\"\"u(m/s)\"\"Press(kPa)\"\"Mach\"\"U1\"\"U2\"\"U3\"\n");
  fprintf(fp3 , "I=%i\n",N);
  fprintf(fp3, "DATAPACKING=POINT\n");
  fprintf(fp3, "DT = (DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )\n");
  
  fp4 = fopen("data/ss/q1Dhistory.dat","w");
  fprintf(fp4, "TITLE = \"Quasi-1D Residual Supersonic History\"\n");
  fprintf(fp4,"variables=\"Iteration\"\"Time(s)\"\"Res1\"\"Res2\"\"Res3\"\n");
  fprintf(fp4, "DATAPACKING=POINT\n");
  fprintf(fp4, "DT = (INTEGER DOUBLE DOUBLE DOUBLE DOUBLE)\n");

  fp5 = fopen("data/sb/q1Dhistory.dat","w");
  fprintf(fp5, "TITLE = \"Quasi-1D Residual Subsonic History\"\n");
  fprintf(fp5,"variables=\"Iteration\"\"Time(s)\"\"Res1\"\"Res2\"\"Res3\"\n");
  fprintf(fp5, "DATAPACKING=POINT\n");
  fprintf(fp5, "DT = (INTEGER DOUBLE DOUBLE DOUBLE DOUBLE)\n");

}


/***************************** END FILE HANDLING *******************************/

