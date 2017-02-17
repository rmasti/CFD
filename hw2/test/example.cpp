#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "nozzle.hpp"

TEST_CASE( "Factorials are computed", "[factorial]" ) {
  double blah = A_x(0.5);
    REQUIRE( blah >= 20.0 );
}
