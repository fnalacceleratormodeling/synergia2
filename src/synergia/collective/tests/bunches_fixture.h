#ifndef IBUNCHES_FIXTURE_H_
#define IBUNCHES_FIXTURE_H_
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "bunch_fixture.h"
#include "synergia/foundation/physical_constants.h"

//from bunch_fixture.h
// const int charge = pconstants::proton_charge;
// const double mass = pconstants::mp;
// const double real_num = 1.0e11;
// const int total_num = 20;
// const double total_energy = 125.0;

//const double mass = pconstants::mp;
//const double total_energy = 1.25*mass;
//int total_num;
//const double real_num = 2.0e12;
const int local_num = 9;
const int turns = 17;
const double turn_length = 246.8;
const double partial_s = 123.4;



// void
// dummy_populate(Bunch &bunch)
// {
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//         for (int i = 0; i < 6; i += 1) {
//             bunch.get_local_particles()[part][i] = 10.0 * part + (1.0 + part
//                     * part / 1000.0) * i + 1000 * bunch.get_comm().get_rank();
//         }
//         bunch.get_local_particles()[part][Bunch::id] = part + total_num
//                 * bunch.get_comm().get_rank();
//     }
// }

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

void
dummy_populate_gaussian(Bunch &bunch, int offset = 0)
{

   Commxx_sptr comm_sptr(new Commxx);
   unsigned long int seed = 718281828;
   Random_distribution distribution(seed, *comm_sptr);

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

struct Bunches_fixture
{
    Bunches_fixture() :
            num_bunches(3), reference_particle(pconstants::electron_charge,
                    mass, total_energy), bunches()
    {
        BOOST_TEST_MESSAGE("setup fixture");
	Commxx_sptr parent_sptr(new Commxx);
        Commxxs commxxs(generate_subcomms(parent_sptr, num_bunches));
        for (int i = 0; i < num_bunches; ++i) {
            Bunch_sptr bs(
                    new Bunch(reference_particle, local_num, real_num,
                            commxxs.at(i)));
            bunches.push_back(
                    bs);
            dummy_populate(*(bunches.back()));
        }
        Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic("dummy.h5"));
	Diagnosticss diagnosticss;
	diagnosticss.push_back(diagnostics_sptr);
	train_diagnosticss.push_back(diagnosticss);
    }

    ~Bunches_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    const int num_bunches;
    Reference_particle reference_particle;
    Bunches bunches;
    Train_diagnosticss  train_diagnosticss;
};



#endif /* IBUNCH_FIXTURE_H_ */
