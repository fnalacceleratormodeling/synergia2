#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/utils/cereal_files.h"

const std::string name("foo");

TEST_CASE("call_evil_function")
{
    Lattice lattice(name);
    REQUIRE_NOTHROW(lattice.do_not_call_this_function());
}
