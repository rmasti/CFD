/* 
 * Computational Fluid Dynamics
 * This code contains the main script for the final project in which all the cases will be handled here
 * Written By: Robert Masti
 * 04/12/2017
 * This is the main file that will that will use a header file containing to function prototypes, structure declarations, symbolic constants, class declaration, etc.. And will use a seperate cpp file for the function declarations. 
 */
#include "fp.hpp"

int main( int argc, char *argv[] )
{
  // Define Constants with the datastructure "constants" 
  const char * InputFile =  argv[1];
  constants C = loadInputFile(InputFile); // from argv load the input constants for data extraction

  // Now with the constants inputted create the output folder for cleanliness sake
  string Output_Folder = buildCaseFolder(C); // return the name of the folder it created for writing to
  Output_Folder = "./output/" + Output_Folder;

  ////////////////////// READ THE MESH ////////////////////////
  MatrixXd xn; // x coord for interface or nodal points of the cells
  MatrixXd yn; // y coord for interface or nodal points of the cells
  MatrixXd zn; // z coord for nodal z points

  MatrixXd xc; // x coord for cell x values
  MatrixXd yc; // y coord for cell y values
  MatrixXd zc; // z coord for cell z values

  // GET COORDS
  inputMesh(xn, yn, zn, xc, yc, zc, C); // resize and fill coords

  // Create xc and yc that includes ghost cells
  MatrixXd xc_g(xc.rows()+2*C.num_ghost, xc.cols()+2*C.num_ghost);
  MatrixXd yc_g(yc.rows()+2*C.num_ghost, yc.cols()+2*C.num_ghost);
  extrapCopyCoords(xc_g, yc_g, xc, yc, C); // get ghost cell coord

  // Create temperature dist for sake of BC
  MatrixXd T(xc_g.rows(), xc_g.cols()); // K

  ////////////////////// SETUP VARIABLES //////////////////////
  // SETUP AREA's
  //i or x dir areas the columns of the matrix
  MatrixXd Ai(xc.rows(), xn.cols()); // sets (row,col) (j,i)
  //j or y dir areas the rows or the matrix
  MatrixXd Aj(xn.rows(), xc.cols()); // sets (row,col) (j,i)

  // SETUP NORMAL VECS
  //Ai normal vector has x and y components (so same size as Ai)
  MatrixXd n_i_xhat(Ai.rows(), Ai.cols());
  MatrixXd n_i_yhat(Ai.rows(), Ai.cols());
  //Aj normal vector has x and y componenets
  MatrixXd n_j_xhat(Aj.rows(), Aj.cols());
  MatrixXd n_j_yhat(Aj.rows(), Aj.cols());
  //Setup History tracking vectors
  MatrixXd L2hist(NEQ, C.nmax); // large matrix with all hist
  MatrixXd Iterhist(NEQ, C.nmax); // same for iter err
  // SETUP Max Speed AND dt
  ////calcs max character..(MHD to) (soundspeed or fast magnetosonic speed)
  MatrixXd MaxSpeed(xc.rows(), xc.cols());
  MatrixXd dt(xc.rows(), xc.cols()); //time step for local or global 
  double dt_min; // used for global t stepping only
  // SETUP VOLUMES 
  MatrixXd Volume(xc.rows(), xc.cols());
  // SETUP NEQ MATRIX ARRAY'S WITH POINTERs
  // direction independent
  MatrixXd* V=new MatrixXd[NEQ]; // primitive variables
  MatrixXd* V_MMS=new MatrixXd[NEQ]; // MMS solution
  MatrixXd* U=new MatrixXd[NEQ]; // conservative variables
  MatrixXd* S=new MatrixXd[NEQ]; // source term
  MatrixXd* Res=new MatrixXd[NEQ]; //residual term only domain of interest
  MatrixXd* U_RK=new MatrixXd[NEQ]; //time step for local or global 
  MatrixXd* Error=new MatrixXd[NEQ]; //calcs error over the domain all eqs
  // direction dependent
  MatrixXd* F=new MatrixXd[NEQ]; // F flux 
  MatrixXd* G=new MatrixXd[NEQ]; // G flux 
  MatrixXd* V_L =new MatrixXd[NEQ]; // Left Prim Var states
  MatrixXd* V_R =new MatrixXd[NEQ]; // Right Prim Var States
  MatrixXd* V_T =new MatrixXd[NEQ]; // Top Prim Var states
  MatrixXd* V_B =new MatrixXd[NEQ]; // Bottom Prim Var States
  // set size of the array 
  // add the ghost cell layers for the big mats
  int ncols = xc.cols()+2*C.num_ghost;
  int nrows = xc.rows()+2*C.num_ghost;
  for (int eq = 0 ; eq < NEQ ; eq++)
  {
    ///// Direction Indepndent /////
    // Interior plus ghost cellvals resize
    U[eq].resize(nrows,ncols); // numerical solution
    V[eq].resize(nrows,ncols); // Numerical solution
    V_MMS[eq].resize(nrows,ncols); // MMS solution
    U_RK[eq].resize(nrows,ncols); //Runge kutta intermediate holder
    // Interior cellvals only resize
    Res[eq].resize(xc.rows(),xc.cols()); // add res only to interior ****
    S[eq].resize(xc.rows(),xc.cols()); // interior only
    Error[eq].resize(xc.rows(),xc.cols()); //error
    ///// Direction Depndent /////
    // Interior facevals xdir only resize
    F[eq].resize(xc.rows(),xn.cols()); // only need flux on the interior
    V_L[eq].resize(xc.rows(),xn.cols()); // all faces of the cells of int
    V_R[eq].resize(xc.rows(),xn.cols()); // right states
    // Interior facevals ydir only resize
    G[eq].resize(xn.rows(),xc.cols()); // only need flux on the interior
    V_T[eq].resize(xn.rows(),xc.cols()); // all faces in y dir to get G
    V_B[eq].resize(xn.rows(),xc.cols()); // all faces in y dir
  }

  ////////////////////// INITIALIZE ///////////////////////////
  //______ SET GEOMETRY ______//
  computeArea(Ai, Aj, xn, yn); // take nodal coords extract interface A
  computeNormalVectors(n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat,
      xn, yn, Ai, Aj); // grab the norm vec 4 computational domain
  computeVolume(Volume, xn, yn); 
  outputArray(Output_Folder,"volume", Volume, 0);
  //______SET PRIMITIVE VARIABLES______//
  cout << "Initializing... " << endl;
  initialize(V, xc_g, yc_g, C);
  //______SET BOUNDARY CONDITIONS/SOURCES______//
  computeTemperature(T, V); // compute temp for slip wall condition
  if (C.f_case == 1) // MMS Curvilinear mesh
  {
    solveSourceMMS(S,xc,yc, C); // grab the source term for the MMS case
    solveSolutionMMS(V_MMS, xc_g, yc_g, C); // solution at cell w/g
    setBCMMS(V, V_MMS, C); // Sets the boundary conditions for MMS
    outputArray(Output_Folder, "S1", S[rhoid], 0); // output Sols
    outputArray(Output_Folder, "S2", S[rhouid], 0);
    outputArray(Output_Folder, "S3", S[rhovid], 0);
    outputArray(Output_Folder, "S4", S[rhoetid], 0);
    outputArray(Output_Folder, "rho_MMS", V_MMS[rhoid], 0);
    outputArray(Output_Folder, "u_MMS", V_MMS[rhouid], 0);
    outputArray(Output_Folder, "v_MMS", V_MMS[rhovid], 0);
    outputArray(Output_Folder, "p_MMS", V_MMS[rhoetid], 0);

  }
  else // case 2 and 3 BC enforce and Source Determ
  {
    // zero the source term only need for MMS
    for (int eq = 0; eq < NEQ; eq++)
      S[eq] = MatrixXd::Constant(xc.rows(), xc.cols(), 0.0);
    // set BC's
    setBC(V, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, T, C);
  }
  primToCons(U, V, C);  U_RK = U; // initialize U_RK
  ////////////////// ENTER INITIALIZATION ////////////////
  MUSCL(V_L,V_R,V_B,V_T,V,C);  // get the left Right Top and Bot states
  compute2DFlux(F,G,n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, V_L, V_R, V_B, V_T, C); //Roe or VL in C
  computeRes(Res, F, G, S, Ai, Aj, Volume, C); // compute Residual knowing fluxes, sources, geometry
  // compute maxspeed
  computeMaxSpeed(MaxSpeed, V, C);
  // compute timestep BOTH Locally and Globally 
  dt_min = computeTimeStep(dt, Volume, Ai, Aj, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, MaxSpeed, V, C);
  // compute norms and errors for the initial step
  L2hist.col(0) = computeL2(Res, C); // Get the initial Residual for normalization
  if (C.f_case == 1) // compute error for MMS case only
  {
    computeError(Error, V, V_MMS, C);
    Iterhist.col(0) = computeL2(Error, C);
  }

  ///////////////////// OUTPUT INITIAL ARRAYS ////////////////////////////
  outputArray(Output_Folder, "rho", V[rhoid], 0);
  outputArray(Output_Folder, "u", V[uid], 0);
  outputArray(Output_Folder, "v", V[vid], 0);
  outputArray(Output_Folder, "p", V[pid], 0);

  outputArray(Output_Folder, "xc", xc, 0);
  outputArray(Output_Folder, "yc", yc, 0);
  outputArray(Output_Folder, "xc_g", xc_g, 0);
  outputArray(Output_Folder, "yc_g", yc_g, 0);


  outputArray(Output_Folder, "Res1", Res[rhoid], 0);
  outputArray(Output_Folder, "Res2", Res[rhouid], 0);
  outputArray(Output_Folder, "Res3", Res[rhovid], 0);
  outputArray(Output_Folder, "Res4", Res[rhoetid], 0);

  outputArray(Output_Folder, "Ai",Ai, 0);
  outputArray(Output_Folder, "Aj",Aj, 0);

  outputArray(Output_Folder, "n_i_x", n_i_xhat, 0);
  outputArray(Output_Folder, "n_i_y", n_i_yhat, 0);
  outputArray(Output_Folder, "n_j_x", n_j_xhat, 0);
  outputArray(Output_Folder, "n_j_y", n_j_yhat, 0);


  outputArray(Output_Folder, "F1", F[frhouid], 0);
  outputArray(Output_Folder, "F2", F[frhouuid], 0);
  outputArray(Output_Folder, "F3", F[frhouvid], 0);
  outputArray(Output_Folder, "F4", F[frhouhtid], 0);

  outputArray(Output_Folder, "G1", G[grhovid], 0);
  outputArray(Output_Folder, "G2", G[grhouvid], 0);
  outputArray(Output_Folder, "G3", G[grhovvid], 0);
  outputArray(Output_Folder, "G4", G[grhovhtid], 0);

  outputArray(Output_Folder, "rho_L", V_L[rhoid], 0);
  outputArray(Output_Folder, "u_L", V_L[uid], 0);
  outputArray(Output_Folder, "v_L", V_L[vid], 0);
  outputArray(Output_Folder, "p_L", V_L[pid], 0);

  outputArray(Output_Folder, "rho_R", V_R[rhoid], 0);
  outputArray(Output_Folder, "u_R", V_R[uid], 0);
  outputArray(Output_Folder, "v_R", V_R[vid], 0);
  outputArray(Output_Folder, "p_R", V_R[pid], 0);

  outputArray(Output_Folder, "rho_B", V_B[rhoid], 0);
  outputArray(Output_Folder, "u_B", V_B[uid], 0);
  outputArray(Output_Folder, "v_B", V_B[vid], 0);
  outputArray(Output_Folder, "p_B", V_B[pid], 0);

  outputArray(Output_Folder, "rho_T", V_T[rhoid], 0);
  outputArray(Output_Folder, "u_T", V_T[uid], 0);
  outputArray(Output_Folder, "v_T", V_T[vid], 0);
  outputArray(Output_Folder, "p_T", V_T[pid], 0);

  /////////////// RUNGE KUTTA ITERATION ///////////////////////

  cout << "Entering Main Time Loop...." << endl;
  for ( int n = 1; n < C.nmax-1; n++) 
  {
    for (int k = 0; k < C.rk_order; k++) 
    {
      // RK takes k steps for a given dt 
      // RK will refine the timestep guess of the residual thru sub cycling
      rungeKutta(U_RK, U, Res, Volume, dt, k, dt_min, C); // update interior
      consToPrim(V, U_RK, C); // convert to prim var and apply BC's
      computeTemperature(T, V); // grab temp for B.C.'s
      if (C.f_case != 1) {
        ////// NEED TO DO BC'S NACA and Inlet
        setBC(V, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, T, C);
      }
      MUSCL(V_L,V_R,V_B,V_T,V,C); // Use new V for extrap to get primvars LRBT
      compute2DFlux(F,G,n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, V_L, V_R, V_B, V_T, C); // get F, G
      computeRes(Res, F, G, S, Ai, Aj, Volume, C);  // Get refined residual and get U_RK again
    }
    // compute maxspeed to get new dt
    computeMaxSpeed(MaxSpeed, V, C);
    // compute timestep local and global with the new vars
    dt_min = computeTimeStep(dt, Volume, Ai, Aj, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, MaxSpeed, V, C);
    U = U_RK; // update RK

    L2hist.col(n+1) = computeL2(Res, C);
    if (C.f_case == 1) // compute error
    {
      computeError(Error, V, V_MMS, C);
      Iterhist.col(n+1) = computeL2(Error, C);
    }

    if ( n % C.pint == 0)
    {
      cout <<  
        "rho: " << L2hist(rhoid,n+1)/L2hist(rhoid,0) << 
        " u: " << L2hist(uid,n+1)/L2hist(uid,0) << 
        " v: " << L2hist(vid,n+1)/L2hist(vid,0) <<  
        " p: " << L2hist(pid,n+1)/L2hist(pid,0) << " n: " << n << endl;
    }

    if ( n % C.wint == 0)
    {
      outputArray(Output_Folder, "rho", V[rhoid], n);
      outputArray(Output_Folder, "u", V[uid], n);
      outputArray(Output_Folder, "v", V[vid], n);
      outputArray(Output_Folder, "p", V[pid], n);

      outputArray(Output_Folder, "Res1", Res[rhoid], n);
      outputArray(Output_Folder, "Res2", Res[rhouid], n);
      outputArray(Output_Folder, "Res3", Res[rhovid], n);
      outputArray(Output_Folder, "Res4", Res[rhoetid], n);
    }

    //if (n==C.nmax/200)
      //C.f_limiter=0; // turn off limiter
    // Check for convergence
    if ( L2hist(rhoid,n+1)/L2hist(rhoid,0) < C.tol && L2hist(uid,n+1)/L2hist(uid,0) < C.tol &&  L2hist(vid,n+1)/L2hist(vid,0) < C.tol && L2hist(pid,n+1)/L2hist(pid,0) < C.tol )
    {
      // output final soln
      outputArray(Output_Folder, "rho", V[rhoid], n);
      outputArray(Output_Folder, "u", V[uid], n);
      outputArray(Output_Folder, "v", V[vid], n);
      outputArray(Output_Folder, "p", V[pid], n);
      outputArray(Output_Folder, "rhou", U[rhouid], n);
      outputArray(Output_Folder, "rhov", U[rhovid], n);
      outputArray(Output_Folder, "rhoe", U[rhoetid], n);

      // Resize history to correct size
      L2hist.conservativeResize(NEQ,n);
      outputArray(Output_Folder, "L2Hist", L2hist, n);
      if (C.f_case == 1) // resize iterhist and output only for case 
      {
        Iterhist.conservativeResize(NEQ,n);
        outputArray(Output_Folder, "IterErrHist", Iterhist, n);
      }
      break;
    }
    if (std::isnan(L2hist(0,n+1)/L2hist(0,0))) // check for NAN's
    {
      cerr << "NAN: 404 abort" << endl;
      exit(1);
    }
  }
  return 0;
}


