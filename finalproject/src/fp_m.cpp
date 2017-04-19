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
  constants C = loadInputFile(InputFile);
  
  // Now with the constants inputted create the output folder for cleanliness sake
 
  string Ouput_Folder = buildCaseFolder(C); // return the name of the folder it created for writing to

  /////////////////////////////////////////////////////////////
  ////////////////////// READ THE MESH ////////////////////////
  /////////////////////////////////////////////////////////////
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
  extrapCopyCoords(xc_g, yc_g, xc, yc, C);
 
  // Create temperature dist for sake of BC
  MatrixXd T(xc_g.rows(), xc_g.cols()); // K
  /////////////////////////////////////////////////////////////
  ////////////////////// SETUP VARIABLES //////////////////////
  /////////////////////////////////////////////////////////////
 
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

  // SETUP Max Speed AND DT
  ////calcs max character..(MHD to) (soundspeed or fast magnetosonic speed)
  MatrixXd MaxSpeed(xc.rows(), xc.cols());

  MatrixXd Dt(xc.rows(), xc.cols()); //time step for local or global 

  // SETUP VOLUMES
  MatrixXd Volume(xc.rows(), xc.cols());

  ////////// matrix arrays
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
  // add the ghost cell layers
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
    V_L[eq].resize(xc.rows(),xn.cols()); // all faces of the cells 
    V_R[eq].resize(xc.rows(),xn.cols());
    // Interior facevals ydir only resize
    G[eq].resize(xn.rows(),xc.cols()); // only need flux on the interior
    V_T[eq].resize(xn.rows(),xc.cols()); // all faces in y dir
    V_B[eq].resize(xn.rows(),xc.cols()); // all faces in y dir
  }

  /////////////////////////////////////////////////////////////
  ////////////////////// INITIALIZE ///////////////////////////
  /////////////////////////////////////////////////////////////

  //______ SET GEOMETRY ______//

  computeArea(Ai, Aj, xn, yn);
  computeNormalVectors(n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat,
      xn, yn, Ai, Aj);
  computeVolume(Volume, xn, yn);


  //______SET PRIMITIVE VARIABLES______//
  cout << "Initializing... " << endl;
  initialize(V, xc_g, yc_g, C);

  computeTemperature(T, V);

  //______SET BOUNDARY CONDITIONS/SOURCES______//
  cout << "Applying BC's... " << endl;
  if (C.f_case == 1) // MMS Curvilinear mesh
  {
    solveSourceMMS(S,xc,yc, C);
    solveSolutionMMS(V_MMS, xc_g, yc_g, C); // solution at cell w/ g
    setBCMMS(V, V_MMS, C);
  }
  else // case 2 and 3 BC enforce and Source Determ
  {
    // zero the source term
    for (int eq = 0; eq < NEQ; eq++)
      S[eq] = MatrixXd::Constant(xc.rows(), xc.cols(), 0.0);
    // set BC STILL TO DO
  }
  /////////////////////////////////////////////////////////////
  //////////////// EXTRAPOLATE, COMPUTE FLUXES ////////////////
  /////////////////////////////////////////////////////////////
  cout << "Performing MUSCL Extrap... " << endl; 
  MUSCL(V_L,V_R,V_B,V_T,V,C); // MATCHES for case 1
  cout << "Computing Fluxes... " << endl; 
  compute2DFlux(F,G,n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, V_L, V_R, V_B, V_T, C);

//  cout << V_L[rhoid].rows() << " c=" << V_L[rhoid].cols() <<  endl;
  /////////////////////////////////////////////////////////////
  /////////////////// OUTPUT ARRAYS ///////////////////////////
  /////////////////////////////////////////////////////////////

  outputArray("n_i_x", n_i_xhat, 0);
  outputArray("n_i_y", n_i_yhat, 0);
  outputArray("n_j_x", n_j_xhat, 0);
  outputArray("n_j_y", n_j_yhat, 0);
 
  
  outputArray("F1", F[frhouid], 0);
  outputArray("F2", F[frhouuid], 0);
  outputArray("F3", F[frhouvid], 0);
  outputArray("F4", F[frhouhtid], 0);
 
  outputArray("G1", G[grhovid], 0);
  outputArray("G2", G[grhouvid], 0);
  outputArray("G3", G[grhovvid], 0);
  outputArray("G4", G[grhovhtid], 0);
  
  outputArray("rho", V[rhoid], 0);
  outputArray("u", V[uid], 0);
  outputArray("v", V[vid], 0);
  outputArray("p", V[pid], 0);

  outputArray("rho_L", V_L[rhoid], 0);
  outputArray("u_L", V_L[uid], 0);
  outputArray("v_L", V_L[vid], 0);
  outputArray("p_L", V_L[pid], 0);
 
  outputArray("rho_R", V_R[rhoid], 0);
  outputArray("u_R", V_R[uid], 0);
  outputArray("v_R", V_R[vid], 0);
  outputArray("p_R", V_R[pid], 0);

  outputArray("rho_B", V_B[rhoid], 0);
  outputArray("u_B", V_B[uid], 0);
  outputArray("v_B", V_B[vid], 0);
  outputArray("p_B", V_B[pid], 0);
 
  outputArray("rho_T", V_T[rhoid], 0);
  outputArray("u_T", V_T[uid], 0);
  outputArray("v_T", V_T[vid], 0);
  outputArray("p_T", V_T[pid], 0);

  outputArray("rho_mms", V_MMS[rhoid], 0);
  outputArray("u_mms", V_MMS[uid], 0);
  outputArray("v_mms", V_MMS[vid], 0);
  outputArray("p_mms", V_MMS[pid], 0);
 
  outputArray("S1", S[rhoid], 0);
  outputArray("S2", S[rhouid], 0);
  outputArray("S3", S[rhovid], 0);
  outputArray("S4", S[rhoetid], 0);


  /////////////////////////////////////////////////////////////
  ////////////////////// COMPUTE FLUX /////////////////////////
  /////////////////////////////////////////////////////////////

  //cout << x << endl;// Identify that all simulations have finished
  return 0;
}


