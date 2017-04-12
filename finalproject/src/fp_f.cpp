/*
 * Computational Fluid Dynamics
 * Written ByRobert Masti
 * 1/27/2017
 * This file contains all function declarations that will be called in hw4_m.cpp which is the 
 * main cpp file. This is the function file, which has its function prototypes stored in the 
 * header file.
 */
#include "fp.hpp" //structure templates and func prototypes


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

  string filename = readMeshName(C);

  cout << " Hellow World 2 " << endl;

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

  system(exec1); // run the bash command with exec
  system(exec2);

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
  return C;
}


