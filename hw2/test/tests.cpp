#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "nozzle.hpp"
#define reltol 1.0e-12
#include <cmath>
using namespace std;



TEST_CASE( "Test deriv Area function at x = 0.11135", "[dAdx]" ) {
  double dA = dAdx(0.11135);
  double predicted = 0.43068128665083;
  REQUIRE( abs(dA - predicted)/predicted <= reltol );
}

TEST_CASE( "Test Area function at x = 0.4215", "[A]" ) {
  double A = A_x(0.4215);
  double predicted = 0.50235087913268;
  REQUIRE( abs(A - predicted)/predicted <= reltol );
}


TEST_CASE( "Test M function at x = -0.21315", "[M]" ) {
  double M = M_x(-0.21315);
  double predicted = 0.808165;
  REQUIRE( abs(M - predicted)/predicted <= reltol );
}
TEST_CASE( "Test timestep calculator", "[dt]" ) {
  constants C;
  C.gamma = 1.4;
  vector<primvar> Vold;
  primvar Vtemp;
  Vtemp.rho = 1.0889;
  Vtemp.u = 435;
  Vtemp.p = 1.55e5; //Pa
  Vold.push_back(Vtemp);
  int i = 0;
  double dt = compute_timestep(Vold, i, C);
  double predicted = 2.2690864001413e-6;
  REQUIRE( abs(dt-predicted)/predicted <= reltol );
}

TEST_CASE( "Volume calculator by avg A's", "[V]" ) {
  vector<double> Xarr;
  Xarr.push_back(-0.3245);
  Xarr.push_back(-0.1234);

  int i = 0;
  vector<double> ALR;
  double volume = compute_volume(Xarr, i, ALR);
  double predicted = 0.0623559347;
  REQUIRE( abs(volume-predicted)/predicted <= reltol*1000 );
}

TEST_CASE( "prim to cons", "[VtoU]" ) {
  primvar Varr;
  Varr.rho = 1.243212312;
  Varr.u = 324.53432;
  Varr.p = 2.92e5;
  constants C;
  C.gamma = 1.4;
  int i = 0;
  consvar U;
  U = primtocons(Varr, C);
  double pred1 = abs(U.rho - Varr.rho)/Varr.rho;
  double pred2 = abs(U.rhou - (Varr.rho*Varr.u))/Varr.rho*Varr.u;
  double temp = (Varr.p)/(C.gamma-1.0);
  temp += 0.5*Varr.rho*Varr.u*Varr.u;
  double pred3 = abs(U.rhoet - temp)/temp;
  double pred = max(pred1,pred2);
  pred = max(pred, pred3);
  //cout << " pred:  " << pred1 << " " << pred2 << " " << pred3 << endl;
  REQUIRE( pred <= reltol );
}

TEST_CASE( "M to prims", "[MtoV]" ) {

  double M = 0.765438987;
  constants C;
  C.gamma = 1.4;
  C.T0 = 2345;
  C.p0 = 3.1e2;
  primvar V = Mtoprim(M, C);
  double psi = 1+0.5*(C.gamma-1)*(M*M);
  double P = (C.p0*1000)*pow(psi, (C.gamma-1)/C.gamma);
  double T = C.T0/psi;
  double rho = P/(R*T);
  double u = sqrt(P*C.gamma/rho)*M;
  double pred1 = abs(V.rho - (rho))/(rho);
  double pred2 = abs(V.u - u)/(u);
  double pred3 = abs(V.p - P)/P;
  double pred = max(pred1,pred2);

  pred = max(pred, pred3);
  //cout << " pred:  " << pred1 << " " << pred2 << " " << pred3 << endl;
  REQUIRE( pred <= reltol );
}
