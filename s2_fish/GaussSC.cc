#include <stdio.h>
#include "macro_bunch_store.h"
#include "basic_toolkit/PhysicsConstants.h"

#include "BasErs_field.h"

using namespace boost::python;

int
apply_BasErs_kick(Macro_bunch_store &mbs, double sigmaX, double sigmaY, double tau)

{

    const double epsilon_0 = PH_MKS_eps0;

    double sigma[2];
    sigma[0] = sigmaX;
    sigma[1] = sigmaY;
    BasErs_field *myfield = new BasErs_field(sigma);

    for (int ipart = 0; ipart < mbs.local_num; ipart++) {
        // Get the field at the particle (x,y) location
        std::vector<double> Efield(3);

        Efield = myfield->NormalizedEField(mbs.local_particles(0, ipart), mbs.local_particles(2, ipart));

        //std::cout << " Ex = " << Efield[0] << " Ey = "<< Efield[1] << " Ez = " << Efield[2] << std::endl;

        // Now x and y kick
        for (int n_axis = 0; n_axis < 2; n_axis++) {
            double gamma = -1 * mbs.ref_particle(5);
            double beta = sqrt(gamma * gamma - 1.0) / gamma;
            const  double c = PH_MKS_c;
            const  double pi = M_PI;

            double mass = mbs.mass * 1.0e9;
            double eps0 = PH_MKS_eps0;
            double Brho = gamma * beta * mass / c;
            double perveance0 = mbs.total_current / (2 * pi * eps0 * Brho * gamma * gamma* \
                                beta * c * beta * c);
            double LengthScale = mbs.units(2 * n_axis);

            double factor = pi * perveance0 * gamma * beta * beta / LengthScale * 4.0 * pi;
            // (We think) this is for the Lorentz transformation of the transverse
            // E field.
            //factor *= gamma;
            // the above fudge seemed to work for high energy, for low energy, I had to multiply current by
            // pi/betagamma (checked using channel).  Fudge factor recalculated below:

            factor *= pi / beta;

            int index = 2 * n_axis + 1; // for n_axis = (0,1,2) Cartesian coordinate x,y,z,
            // in particle store indexing, px,py,pz = (1,3,5)
            double kick = Efield[n_axis] * tau * factor ; // the SC kick with the 2D geometric field contribution

            //if (ipart < 5) {
            //std::cout << " our kick " << kick << "the 3d kick " << kick2 << " ID " << mbs.local_particles(6,ipart) << std::endl;}

            // update the data structure
            mbs.local_particles(index, ipart) += kick;
        }

    }

    delete myfield;

    return 0;
}

BOOST_PYTHON_MODULE(GaussSC)
{
    def("apply_BasErs_kick", &apply_BasErs_kick);
}
