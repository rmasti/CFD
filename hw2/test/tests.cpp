#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "nozzle.hpp"
#define reltol 1.0e-12
using namespace std;
TEST_CASE( "Test Area function at x = 0.4215", "[A]" ) {
  double A = A_x(0.4215);
  double predicted = 0.50235087913268;
  REQUIRE( abs(A - predicted)/predicted <= reltol );
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
  double predicted = 1.1345432000707e-6;
  REQUIRE( abs(dt-predicted)/predicted <= reltol );
}

TEST_CASE( "Volume calculator by avg A's", "[V]" ) {
  vector<double> Xarr;
  Xarr.push_back(-0.3245);
  Xarr.push_back(-0.1234);

  int i = 0;
  
  double volume = compute_volume(Xarr, i);
  double predicted = 0.0623559347;
  cout << " pred " << predicted << " act " << volume << endl;
  REQUIRE( abs(volume-predicted)/predicted <= reltol*1000 );
}
