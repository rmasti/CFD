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

  inputMesh(xn, yn, zn, xc, yc, zc, C); // resize and fill coords

  /////////////////////////////////////////////////////////////
  ////////////////////// INITIALIZE ///////////////////////////
  /////////////////////////////////////////////////////////////


  //cout << x << endl;// Identify that all simulations have finished
  return 0;
}


