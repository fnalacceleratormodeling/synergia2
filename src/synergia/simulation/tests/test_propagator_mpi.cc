#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
#include "synergia/bunch/diagnostics_particles.h"


BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

struct Object_to_sptr_hack
{
    void
    operator()(void const *) const
    {
    }
};



 BOOST_FIXTURE_TEST_CASE(propagate, Lattice_fixture)
 {
   Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
           "space_charge"));

   Commxx mcommx=Commxx();
   Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Reference_particle reference_particle=b.bunch.get_reference_particle() ;
    /// bunches with different comms has to be generated
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
    ///*********************************
//    Bunch_with_diagnostics_train bunch_diag_train(num_bunches, bunch_separation, mcommx);
//    for (int bunchnum = 0; bunchnum < num_bunches; ++bunchnum) {
//         if (bunch_diag_train.is_on_this_rank(bunchnum)){
//            Commxx commx=bunch_diag_train.get_comm(bunchnum);
//            Bunch_sptr bunch_sptr(new Bunch(reference_particle, total_num,real_num, commx, particle_charge));
//            Random_distribution distribution(0,commx);
//            populate_6d(distribution, *bunch_sptr, means, covariances);
//
//            Diagnostics_actions_sptr diagnostics_actions_sptr(new Diagnostics_actions);
//            Bunch_with_diagnostics_sptr bunch_with_diagnostics_sptr(new Bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr));
//
//            std::string ss =boost::lexical_cast<std::string>(bunchnum);
//            Diagnostics_sptr step_full2_sptr(new  Diagnostics_full2(bunch_sptr,"full2_per_step_"+ss+".h5"));
//            Diagnostics_sptr turn_particles_sptr(new  Diagnostics_particles(bunch_sptr,"particles_per_turn_"+ss+".h5"));
//            bunch_with_diagnostics_sptr->add_per_step_diagnostics(step_full2_sptr);
//            bunch_with_diagnostics_sptr->add_per_turn_diagnostics(turn_particles_sptr);
//            bunch_diag_train.set_bunch_diag_sptr(bunchnum, bunch_with_diagnostics_sptr);
//         }
//    }
//
//    int num_turns = 4;
//    propagator.propagate(bunch_diag_train, num_turns);

 }

