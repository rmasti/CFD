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
  ////////////////////// SET GEOMETRY /////////////////////////
  /////////////////////////////////////////////////////////////
  MatrixXd xn; // x coord for interface or nodal points of the cells
  MatrixXd yn; // y coord for interface or nodal points of the cells
  MatrixXd zn; // z coord for nodal z points
 
  MatrixXd xc; // x coord for cell x values
  MatrixXd yc; // y coord for cell y values
  MatrixXd zc; // z coord for cell z values

  // GET COORDS
  inputMesh(xn, yn, zn, xc, yc, zc, C); // resize and fill coords

  // SETUP AREA's
  //i or x dir areas the columns of the matrix
  MatrixXd Ai(xc.rows(), xn.cols()); // sets (row,col) (j,i)
  //j or y dir areas the rows or the matrix
  MatrixXd Aj(xn.rows(), xc.cols()); // sets (row,col) (j,i)

  // SETUP NORMAL VECS
  //Ai normal vector has x and y components
  MatrixXd n_i_xhat(xc.rows(), xn.cols());
  MatrixXd n_i_yhat(xc.rows(), xn.cols());
  //Aj normal vector has x and y componenets
  MatrixXd n_j_xhat(xn.rows(), xc.cols());
  MatrixXd n_j_yhat(xn.rows(), xc.cols());

  // SETUP VOLUMES
  MatrixXd Volumn(xc.rows(), xc.cols());

  /////////////////////////////////////////////////////////////
  ////////////////////// INITIALIZE ///////////////////////////
  /////////////////////////////////////////////////////////////
 
  // CONSTRUCT VARS
  // direction independent
  MatrixXd* V=new MatrixXd[NEQ]; // primitive variables
  MatrixXd* U=new MatrixXd[NEQ]; // conservative variables
  MatrixXd* S=new MatrixXd[NEQ]; // source term
  MatrixXd* Res=new MatrixXd[NEQ]; //residual term only on domain of interest
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
    // Interior cellvals only resize
    Res[eq].resize(xc.rows(),xc.cols()); // add res only to interior ****
    S[eq].resize(xc.rows(),xc.cols()); // interior only
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

  cout << " Initialized " << endl;

  //cout << x << endl;// Identify that all simulations have finished
  return 0;
}


