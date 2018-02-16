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

    double c0[4];

    c0[0] = 0.001;
    c0[2] = -0.002;
    c0[1] = 0.0;
    c0[3] = 0.0;

    // c++
    NonlinearLensPropagatorCmplx(knll, cnll, c0);

    for (int i=0; i<4; ++i) {
        if (i != 0) std::cout << " ";
        std::cout << std::setprecision(18) << c0[i];
    }
    std::cout << std::endl;
    return 0;
}

        
