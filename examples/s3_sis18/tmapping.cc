
#include "synergia/foundation/trigon.h"
#include <iostream>


void evaluation()
{
    using Tri = Trigon<double, 3, 3>;

    // coordinates
    Tri x(0.0, 0);
    Tri y(0.0, 1);
    Tri z(0.0, 2);

    std::cout << "x = " << x << "\n";
    std::cout << "y = " << y << "\n";

    // functions
    auto f = exp(x)*exp(y);
    auto g = kt::qpow(f, 3);
    auto k = x + y + x*y + y*z + z*z;

    std::cout << "f = exp(x)*exp(y) = " << f << "\n";
    std::cout << "g = qpow(f,3) = " << g << "\n";
    std::cout << "k = x+y+xy+yz+zz = " << k << "\n";

    // evaluate functions at this coordinate 
    arr_t<double, 3> v{0.1, 0.2, 1.0};

    std::cout << "f(v) = " << f(v) << "\n";
    std::cout << "g(v) = " << g(v) << "\n";
    std::cout << "k(v) = " << k(v) << "\n\n";
}

void composition()
{
    using Tri = Trigon<double, 3/*power*/, 2/*dim*/>;

    // coordinates mapping a
    Tri x(0.0, 0);
    Tri y(0.0, 1);

    std::cout << "x = " << x << "\n";
    std::cout << "y = " << y << "\n";

    // mapping of a
    TMapping<Tri> a;
    a[0] = x*y*y + exp(x+y);
    a[1] = cos(y*x*x) / (x+2.0);

    std::cout << "a[0] = " << a[0] << "\n";
    std::cout << "a[1] = " << a[1] << "\n";

    // coordinates for mapping b
    Tri xx(a[0].value(), 0);
    Tri yy(a[1].value(), 1);

    // mapping b
    TMapping<Tri> b;
    b[0] = sin(xx) * cos(yy);
    b[1] = exp(xx*xx*xx) / (xx*yy);

    std::cout << "b[0] = " << b[0] << "\n";
    std::cout << "b[1] = " << b[1] << "\n";

    // compose b(a) with a reference point
    auto c = b(a, {xx.value(), yy.value()});

    std::cout << "c[0] = " << c[0] << "\n";
    std::cout << "c[1] = " << c[1] << "\n";

    // direct concat
    Tri q = x*y*y + exp(x+y);
    Tri v = cos(y*x*x) / (x+2.0);;

    Tri w = sin(q)*cos(v);
    Tri z = exp(q*q*q) / (q*v);

    std::cout << "w = " << w << "\n";
    std::cout << "z = " << z << "\n";

    return;
}

int main(int argc, char ** argv)
{
    //MPI_Init(&argc, &argv);
    //Kokkos::initialize(argc, argv);

    evaluation();
    composition();

    //Kokkos::finalize();
    //MPI_Finalize();
    return 0;
}

