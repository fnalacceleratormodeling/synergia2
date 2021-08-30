
#define BOOST_CHECK_EQUAL(a, b)     CHECK((a)==(b))
#define BOOST_CHECK_CLOSE(a, b, t)  CHECK((a)==Approx((b)).margin(t))
#define BOOST_CHECK_NO_THROW(a)     REQUIRE_NOTHROW((a))
#define BOOST_CHECK_THROW(a, b)     REQUIRE_THROWS_AS((a), b)


