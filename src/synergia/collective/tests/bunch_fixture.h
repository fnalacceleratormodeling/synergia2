#ifndef BUNCH_FIXTURE_H_
#define BUNCH_FIXTURE_H_
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_print.h"

const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.0e11;
const int total_num = 20;
const double total_energy = 125.0;
struct Bunch_fixture
{
    Bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm(), bunch(reference_particle,
                total_num, real_num, comm), distribution(0, comm)
    {
        BOOST_TEST_MESSAGE("setup bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1);
            }
        }
        populate_6d(distribution, bunch, means, covariances);
    }

    ~Bunch_fixture()
    {
        BOOST_TEST_MESSAGE("teardown bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
    Random_distribution distribution;
};

// const int charge = pconstants::proton_charge;
// const double mass = pconstants::mp;
// const double real_num = 1.0e11;
// const int total_num = 20;
// const double total_energy = 125.0;
// 
// 
// BOOST_GLOBAL_FIXTURE(MPI_fixture)



// const int total_num1 = 1000;
// const double real_num1 = 2.0e12;
// const double step_length = 1.23;
// struct  Bunch_fixture_impedance
// {
//     Bunch_fixture_impedance() :
//         four_momentum(mass, total_energy), reference_particle(
//                 pconstants::proton_charge, four_momentum),
//                 comm(MPI_COMM_WORLD),
//                  bunch(reference_particle, total_num1,
//                         real_num1, comm), step(step_length)
//     {
//         BOOST_TEST_MESSAGE("setup fixture");
//     }
//     ~Bunch_fixture_impedance()
//     {
//         BOOST_TEST_MESSAGE("teardown fixture");
//     }
// 
//     Four_momentum four_momentum;
//     Reference_particle reference_particle;
//     Commxx comm;
//     Bunch bunch;
//     Step step;
// };

void
dummy_populate(Bunch &bunch, int offset = 0)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // coordinates
        for (int i = 0; i < 6; i += 2) {
            bunch.get_local_particles()[part][i] = 10.0 * (part + offset) + i;
        }
        bunch.get_local_particles()[part][4] /=100.;
        // momenta
        for (int i = 1; i < 6; i += 2) {
            bunch.get_local_particles()[part][i] = 1e-4 * (10.0 * (part
                    + offset) + i);
        }
    }
}




#endif /* BUNCH_FIXTURE_H_ */
