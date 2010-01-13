#include <iostream>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

double
hello(double x)
{
    std::cout << "hello world from the C++ function hello.\n";
    std::cout << "hello was passed the argument " << x
        << " and will return " << x+1.0 << std::endl;
    return x+1.0;
}

using namespace boost::python;

BOOST_PYTHON_MODULE(bp_hello)
{
    def("hello",&hello);
}

