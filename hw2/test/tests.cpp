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
  primvar Vtemp;
  Vtemp.rho = 1.0889;
  Vtemp.u = 435;
  Vtemp.p = 1.55e5; //Pa
  vector<consvar> Uold;
  Uold.push_back(primtocons(Vtemp, C));
  int i = 0;
  double dt = compute_timestep(Uold, i, C);
  double predicted = cfl*dx/(881.41200787923);
  REQUIRE( abs(dt-predicted)/predicted <= reltol );
}

TEST_CASE( "Volume calculator by avg A's", "[V]" ) {
  vector<double> Xarr;
  Xarr.push_back(-0.3245);
  Xarr.push_back(-0.1234);

  int i = 0;
  vector<double> ALR(2);
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

TEST_CASE( "Extrapolation test", "[extrapolation]" ) {

  constants C;
  vector<consvar> Uarr(9);
  Uarr[0].rho = 98.321;
  Uarr[1].rho = 9.321;
  Uarr[2].rho = 0.321;
  Uarr[3].rho = 5.321;
  Uarr[4].rho = 3.321;
  Uarr[5].rho = 1.321;
  Uarr[6].rho = 9.321;
  Uarr[7].rho = 21.321;
  Uarr[8].rho = 10.321;

  double Uendpredict = -4.679;

  double Ubeginpredict = 11.321;

  extrapolate_to_ghost(Uarr, C);
  double pred1 = abs(Uarr.front().rho-Ubeginpredict)/abs(Ubeginpredict);

  double pred2 = abs(Uarr.back().rho-Uendpredict)/abs(Uendpredict);

  cout << "val " << Uarr.front().rho << " " << Ubeginpredict << endl;
  cout << "PRED PRED PRED " << pred1 << " " << pred2 << endl;
  double predmax = max(pred1, pred2);
  REQUIRE( predmax <= reltol );
}


TEST_CASE( "constoprim", "[UtoV]" ) {

  constants C;
  C.gamma = 1.4;
  consvar U;
  U.rho = 1.21123;
  U.rhou = 321.456;
  U.rhoet = 678392.0;
 
  primvar Vout = constoprim(U, C);

  double predrho = 1.21123;
  double predu = 265.39633265358;
  double predp = 254294.1512981;
 
  double pred1 = abs(predrho-Vout.rho)/predrho;
  double pred2 = abs(predu-Vout.u)/predu;
  double pred3 = abs(predp-Vout.p)/predp;
 
  double predmax = max(pred1, pred2);
  predmax = max(predmax, pred3);
  REQUIRE( predmax <= reltol );
}

TEST_CASE( "fluxcalc", "[flux]" ) {

  constants C;
  C.gamma = 1.4;
  consvar U;
  U.rho = 1.21123;
  U.rhou = 321.456;
  U.rhoet = 678392.0;
  double rho = U.rho;

  fluxes F = fluxcalc(U, C);

  double u = U.rhou/U.rho;

  double p = (C.gamma - 1.0)*(U.rhoet-0.5*U.rho*u*u);
  
  double predrhou = U.rhou;

  double predrhouu_and_p = U.rho*u*u + p;

  //double predrhouht = U.rho*u*((C.gamma/(C.gamma-1.0))*(p/U.rho) - 0.5*u*u); 

  double predrhouht = U.rhou*((C.gamma/(C.gamma-1.0))*p/rho + 0.5*u*u);


  //cout << predrhouht << " " << F.rhouht << endl;

  double pred1 = abs(predrhou-F.rhou)/predrhou;
  double pred2 = abs(predrhouu_and_p-F.rhouu_and_p)/predrhouu_and_p;
  double pred3 = abs(predrhouht-F.rhouht)/predrhouht;

  //cout << pred1 << " " << pred2 << " " << pred3 << endl;
  double predmax = max(pred1, pred2);
  predmax = max(predmax, pred3);
  REQUIRE( predmax <= 0.001 );
}

