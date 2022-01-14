#include "synergia/utils/catch.hpp"
#include "synergia/lattice/mx_expr.h"
#include "synergia/lattice/mx_parse.h"

using namespace synergia;

const double tolerance = 1.0e-12;

TEST_CASE("mx_expr")
{
    {
        mx_expr expr;
        CHECK(parse_expression("3+1", expr));
        REQUIRE_NOTHROW(mx_eval(expr));
        CHECK(mx_eval(expr) == Approx(4.0).margin(tolerance));
        CHECK(mx_expr_is_number(expr));
    }

    {
        mx_expr expr;
        CHECK(parse_expression("3+a", expr));
        REQUIRE_THROWS(mx_eval(expr));
        CHECK(not mx_expr_is_number(expr));
    }

}


