#include <iostream>

#include <boost/python.hpp>

using namespace boost::python;

class Hello
{
public:
    Hello() {}
    void sayit() { std::cout << "hello\n"; }
    int returnit() { return 459045; } // 459045 is integer for "hello"
};

int
hello()
{
    std::cout << "hello\n";
    return 459045;
}

BOOST_PYTHON_MODULE(pyhello)
{
    def("hello", hello);
    class_<Hello>("Hello", init<>())
            .def("sayit", &Hello::sayit)
            .def("returnit", &Hello::returnit)
            ;
}
