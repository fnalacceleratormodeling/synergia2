#ifndef BUNCH_FIXTURE_H_
#define BUNCH_FIXTURE_H_
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics.h"
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
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                total_num, real_num, comm_sptr), distribution(0, *comm_sptr),
                empty_diagnostics()
    {
    BOOST_TEST_MESSAGE("setup bunch fixture");
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] =0.;// i * 0.0000072;
        for (int j = i; j < 6; ++j) {
            covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1)*0.00000000001;
        }
         covariances[i][i] *= 10.0; // this makes for a positive-definite matrix
    }

    for (int i = 0; i < 6; ++i) {
         covariances[1][i] *=0.01;
         covariances[i][1] *=0.01;
         covariances[3][i] *=0.01;
         covariances[i][3] *=0.01;
         covariances[5][i] *=0.01;
         covariances[i][5] *=0.01;
    }




     populate_6d(distribution, bunch, means, covariances);

    }

    ~Bunch_fixture()
    {
        BOOST_TEST_MESSAGE("teardown bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    Random_distribution distribution;
    Diagnosticss empty_diagnostics;
};

#endif /* BUNCH_FIXTURE_H_ */
