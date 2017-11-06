#include "catch.hpp"
#include "fvCFD.H"

#include "dummy.H"

TEST_CASE("hello")
{
    dummy d;
    d.hello();
}
