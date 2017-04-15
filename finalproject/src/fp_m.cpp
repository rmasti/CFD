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
  MatrixXd* U=new MatrixXd[NEQ]; // conservative variables
  MatrixXd* S=new MatrixXd[NEQ]; // source term
  MatrixXd* Res=new MatrixXd[NEQ]; //residual term only domain of interest
  MatrixXd* U_RK=new MatrixXd[NEQ]; //time step for local or global 
  MatrixXd* Error=new MatrixXd[NEQ]; //calcs error over the domain all eqs

  MatrixXd* Psi_Pos=new MatrixXd[NEQ]; // Psi_pos for finding Left and Right state
  MatrixXd* Psi_Neg=new MatrixXd[NEQ]; // Psi_neg for finding left and right state

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
    U[eq].resize(nrows,ncols);
    V[eq].resize(nrows,ncols);
    U_RK[eq].resize(nrows,ncols); //Runge kutta intermediate holder
    // Interior cellvals only resize
    Res[eq].resize(xc.rows(),xc.cols()); // add res only to interior ****
    S[eq].resize(xc.rows(),xc.cols()); // interior only
    Error[eq].resize(xc.rows(),xc.cols()); //error
    // Interior plus ghost facevals resize
    Psi_Neg[eq].resize(nrows+1,ncols+1); //Also can have negative values
    Psi_Pos[eq].resize(nrows+1,ncols+1); //Every interface can have a psi value*****

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
  initialize(V, xc_g, yc_g, C);
  computeTemperature(T, V);
  outputArray("rho", V[rhoid], 0); // write to the debug file

  //cout << " Initialized " << endl;

  //cout << x << endl;// Identify that all simulations have finished
  return 0;
}


