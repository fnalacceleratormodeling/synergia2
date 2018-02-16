#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

//
// read in particle data from particles.txt
//
// compare tracking of particles through nonlinear lens between
// Chad Mitchell's Fortran routine NonlinearLensPropagatorCmplx
// and NonlinearLensPropagatorCmpx in C++.
//

void NonlinearLensPropagatorCmplx(const double&, const double&, double *);

int main(int argc, char *argv[])
{
    const std::string fname("particles.txt");
    const double tolerance = 1.0e-15;

    double coords[4];
    double chefcoords[4];

    // these values come from the element n.10 and n.11 in
    // the file lattice_1IO_nll_center.madx
    const double knll = 5.479576037e-06;
    const double cnll = 0.008105461952;

    std::fstream partin(fname.c_str(), std::fstream::in);
    std::fstream newpartin("newparticles.txt", std::fstream::in);

    int pnum = 0;
    while (partin.good()) {
        partin >> coords[0];
        partin >> coords[1];
        partin >> coords[2];
        partin >> coords[3];
        ++pnum;

        newpartin >> chefcoords[0];
        newpartin >> chefcoords[1];
        newpartin >> chefcoords[2];
        newpartin >> chefcoords[3];

        double c0[4];

        c0[0] = coords[0];
        c0[2] = coords[2];
        c0[1] = 0.0;
        c0[3] = 0.0;

        // c++
        NonlinearLensPropagatorCmplx(knll, cnll, c0);

        // force failure by uncommenting line below
        //c0[1] += 1.0e-15;

        if (std::abs(c0[1] - chefcoords[1]) > tolerance) {
            std::cout << "difference for particle " << pnum << std::endl;
            std::cout << "c0[1] <--> chefcoords[1]: " << std::setprecision(15) << c0[1] << " <--> " << chefcoords[1] << std::endl;
        }
        if (std::abs(c0[3] - chefcoords[3]) > tolerance) {
            std::cout << "difference for particle " << pnum << std::endl;
            std::cout << "c0[3] <--> chefcoords[3]: " << std::setprecision(15) << c0[3] << " <--> " << chefcoords[3] << std::endl;
        }

//std::cout << pnum << std::endl;
    }
    return 0;
}

        
