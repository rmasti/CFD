/*
 * Computational Fluid Dynamics
 * Written ByRobert Masti
 * 1/27/2017
 * This file contains all function declarations that will be called in hw4_m.cpp which is the 
 * main cpp file. This is the function file, which has its function prototypes stored in the 
 * header file.
 */
#include "fp.hpp" //structure templates and func prototypes



void outputArray(
    // This function will output any eigen matrix into a file
    string FileName,                  //input - File
    MatrixXd& out,                    //input - matrix
    int n)                            //input - iteration
{
  ostringstream StrConvert;
  StrConvert << n; // import n as an ostringstream
  string num_iter = StrConvert.str();// convert to string
  string Address = "./debug/";// add together
  string Suffix = ".txt";
  Address = Address + FileName + num_iter + Suffix; // can combine string
  ofstream outfile; // output
  outfile.open(Address.c_str()); // access string component
  //outfile << setprecision(14) << out.transpose() << endl; // output trans
  outfile << setprecision(14) << out << endl; // output trans
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
    cerr << "ERROR: Size Mismatch in Initialization?!?!" << endl;
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
  //zero out the corners for sake of clarity.
  for (int i = 0; i < C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bot left
      xc_g(j,i) = 0.0; yc_g(j,i) = 0.0;
      // bot right
      xc_g(j,right+i) = 0.0; yc_g(j,right+i) = 0.0;
      // top right
      xc_g(top+j,right+i) = 0.0; yc_g(top+j,right+i) = 0.0;
      // top left
      xc_g(top+j,i) = 0.0; yc_g(top+j,i) = 0.0;
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
  return C;
}


