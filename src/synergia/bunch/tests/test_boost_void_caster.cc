#define BOOST_TEST_MAIN
#include "boost/test/unit_test.hpp"
#include <cstddef>

class Base
{
public:
  virtual ~Base() {}; 
  virtual int foo() = 0;

  int i() const { return i_; }
private:
  int i_ = 1;
};

class Derived : public Base
{
  int foo() override;
};

int Derived::foo() { return i() * 2; };

BOOST_AUTO_TEST_CASE(test1)
{
  Base*    pb = reinterpret_cast<Base*>(1 << 20);
 Derived* pd = static_cast<Derived*>(pb);
 auto dist = reinterpret_cast<std::ptrdiff_t>(pd);
 BOOST_CHECK(dist == (1 << 20));
}
