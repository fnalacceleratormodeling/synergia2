/*******************************************
** main.cc
** Contains:
**
*******************************************/
#include "scalar_field.h"
#include "TripleT.h"
#include <vector>
int main( )
{

    Real_scalar_field my_Scalar_Field( );
    int3 t(1, 2, 3);
    std::cout << "t[0] = " << t[0] << std::endl;
    t[0] = 5;
    std::cout << "new t[0] = " << t[0] << std::endl;
    int foo[3];
    foo[0] = 10;
    foo[1] = 20;
    foo[2] = 30;
    int3 u(foo);
    std::cout << "u[2] = " << u[2] << std::endl;
    std::vector<int> bar(3);
    bar[0] = 100;
    bar[1] = 200;
    bar[2] = 300;
    int3 v(bar);
    std::cout << "v[2] = " << v[2] << std::endl;


    //my_grid.do_something();

    return 0;
}
