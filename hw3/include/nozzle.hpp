/*
 * Comp. Fluid Dynamcis
 * Written By: Robert Masti
 * 1/27/2017
 * This is a header file for nozzlem.cpp which will contain calls from the nozzlef.cpp file which will contain universal constants as well as function prototypes, and even structure declarations.
 */
#ifndef nozzle_H_
#define nozzle_H_
#define N 10
#define xmax 1.5
#define xmin -1.5
#define xmax_dom 1.0
#define xmin_dom -1.0
#define num_ghost_cells 3
#define R 287.058
#define kappa2 (1.0/3.0)
#define kappa4 (1.0/48.0)
#define dx (xmax-xmin)/(N)
#define localdt false
#define nmax 10000000
#define mymax(a,b) ((a>b)?a:b)
#define mymin(a,b) ((a<b)?a:b)

/////////////////////////////////////////////////////////////////////////
///////////////////////// STRUCTURE DEFINITIONS /////////////////////////
/////////////////////////////////////////////////////////////////////////

// prim_exact is a structure that containing structure types double 
struct prim_exact
{
  double rho;
  double u;
  double p;
  double M;
};

struct primvar
{
  double rho;
  double u;
  double v;
  double p;
};

// constants contains doubles and booleans
struct constants
{
  double p0;
  double T0;
  double A_t;
  bool cond; //0 for subsonic vs //1 for 
  double tol;
  double gamma;
  bool outflow;
  double pb; //back pressure
  double cfl;
};
// create consvar structure uses type double
struct consvar
{
  double rho;
  double rhou;
  double rhov;
  double rhoet;
};
struct fluxes
{
  double rhou;
  double rhouu_and_p;
  double rhouht;
};

/////////////////////////////////////////////////////////////////////////
///////////////////////// FUNCTIONS PROTOTYPES //////////////////////////
/////////////////////////////////////////////////////////////////////////
primvar constoprim(
    consvar U, 
    constants const C
    );

consvar primtocons(
    primvar V, 
    constants const C
    );

prim_exact exactsol(
    double A_x, 
    constants const C
    );

primvar Mtoprim(
    double M, 
    constants const C
    );

void output_file_headers();
void quasi1Dnozzle(constants C);
void isentropic_exact(constants C);


void set_geometry(
    std::vector<double> &Aarr, 
    std::vector<double> & xinterface, 
    std::vector<double> & xcenter, 
    std::vector<double> &Marr
    );

double A_x(double xcoord);
double dAdx(double x);
double M_x(double x);




void initialize(std::vector<consvar> &U, std::vector<double> const &M, constants const C);


void write_out(FILE* &fp2, std::vector<double> const &Aarr, std::vector<double> const &XCarr, std::vector<consvar> const &U, constants const C);

void write_res(FILE* &file, std::vector<int> const &n, std::vector<double> const &t, std::vector<double> const &Res1, std::vector<double> const &Res2, std::vector<double> const &Res3);


void set_boundary_cond(std::vector<consvar> &U, constants const C);

void compute_fluxes(std::vector<fluxes> &F, std::vector<consvar> const &U, constants const C);

void reconstruct_U(std::vector<consvar> &U_avg, std::vector<consvar> const &U);

fluxes fluxcalc( consvar const &U, constants const C);

void iteration_step(
    std::vector<consvar> &Unew,         //Output (fill) Conserved Variable std::vector
    std::vector<consvar> &Resarr,       //Output (fill) residual array
    double &dt,
    std::vector<double> &Marr,          //Output (fill) Mach array
    std::vector<fluxes> const &F,       //Input fluxes (with artificial diss)
    std::vector<consvar> const  &Uold,  //Input the old consvar std::vector
    std::vector<double> const &XCarr,   //Input the center cell coords
    std::vector<double> const &Xarr,    //Input the interface coords
    constants const C              //Input constants for conversions, etc.
    );



double compute_timestep(std::vector<consvar> const &Uold, int i, constants const C);

double compute_volume(std::vector<double> &ALR, std::vector<double> const &Xarr, int i );

double primtoM(primvar V, constants const C);


void compute_norms(
    std::vector<double> &Linfnorm,
    std::vector<double> &L1norm,
    std::vector<double> &L2norm,
    std::vector<consvar> const &Resarr
    );



void artificial_viscosity(std::vector<fluxes> &dvec, std::vector<consvar> const &U, constants const C);

fluxes compute_dflux( double const &epsilon2, double const &epsilon4, double const &lambda, int const   &i, std::vector<consvar> const &U);

void extrapolate_to_ghost(std::vector<consvar> &Uarr, constants const C);

#endif

