#ifndef PROPAGATOR_FIXTURE_H_
#define PROPAGATOR_FIXTURE_H_

#include "lattice_fixture.h"
namespace kludge
{
#include "synergia/bunch/tests/bunches_fixture.h"
}

struct Propagator_fixture
{
    Propagator_fixture() :
            l(), bs(), space_charge(
                    new Dummy_collective_operator("space_charge")), lattice_simulator(
                    l.lattice_sptr, 2), stepper_sptr(
                    new Split_operator_stepper(lattice_simulator, space_charge,
                            7)), propagator(stepper_sptr), covariances(
                    boost::extents[6][6]), means(boost::extents[6]), seed(67), distribution(
                    seed, *l.b.comm_sptr)

    {
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.; // i * 0.0000072;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1)
                        * 0.00000000001;
            }
            covariances[i][i] *= 10.0; // this makes for a positive-definite matrix
        }

        for (int i = 0; i < 6; ++i) {
            covariances[1][i] *= 0.01;
            covariances[i][1] *= 0.01;
            covariances[3][i] *= 0.01;
            covariances[i][3] *= 0.01;
            covariances[5][i] *= 0.01;
            covariances[i][5] *= 0.01;
        }
    }
    ;

    ~Propagator_fixture()
    {
        BOOST_TEST_MESSAGE("teardown propagator fixture");
    }
    ;

    Lattice_fixture l;
    kludge::Bunches_fixture bs;
    Dummy_collective_operator_sptr space_charge;
    Lattice_simulator lattice_simulator;
    Split_operator_stepper_sptr stepper_sptr;
    Propagator propagator;
    MArray2d covariances;
    MArray1d means;
    int seed;
    Random_distribution distribution;

};

#endif /* PROPAGATOR_FIXTURE_H_ */
