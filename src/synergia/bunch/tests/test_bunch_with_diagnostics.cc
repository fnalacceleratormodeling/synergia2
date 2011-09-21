#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/bunch_with_diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/tests/bunch_fixture.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture)

// const double tolerance = 1.0e-14;
// 
// const double mass = 100.0;
// const double total_energy = 125.0;
// const int total_num = 100;
// const double real_num = 2.0e12;


BOOST_FIXTURE_TEST_CASE(construct,Bunch_fixture)
{    
    Bunch_sptr bunch_sptr(new Bunch(bunch));

    
     Multi_diagnostics full2_diagnostics;
     Multi_diagnostics particle_diagnostics;
 
     Diagnostics_sptr step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"full2_per_step.h5"));
     full2_diagnostics.append(step_full2_diag_sptr);
    
    Diagnostics_sptr second_step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"second_full2_per_step.h5"));
    full2_diagnostics.append(second_step_full2_diag_sptr);
     
     Diagnostics_sptr turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"particles_per_turn.h5"));
     particle_diagnostics.append(turn_particles_diag_sptr);
   
     Diagnostics_sptr second_turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"second_particles_per_turn.h5"));
     particle_diagnostics.append(second_turn_particles_diag_sptr);
    
     
     Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr,full2_diagnostics, particle_diagnostics);
}

BOOST_FIXTURE_TEST_CASE(construct1,Bunch_fixture)
{    
    Bunch_sptr bunch_sptr(new Bunch(bunch));
    Diagnostics_sptr step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"test1_full2_per_step.h5"));
    Diagnostics_sptr turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"test1_particles_per_turn.h5"));

    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr,step_full2_diag_sptr, turn_particles_diag_sptr);
}

BOOST_FIXTURE_TEST_CASE(construct2,Bunch_fixture)
{    
    Bunch_sptr bunch_sptr(new Bunch(bunch));
    Diagnostics_sptr step_diag_sptr(new  Diagnostics_basic(bunch_sptr,"test_per_step.h5"));
    Diagnostics_sptr turn_diag_sptr(new  Diagnostics_basic(bunch_sptr,"test_per_turn.h5"));

    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr,step_diag_sptr, turn_diag_sptr);
}

/*BOOST_FIXTURE_TEST_CASE(construct2,Bunch_fixture)
{    
   

    
     Multi_diagnostics full2_diagnostics;
     Multi_diagnostics particle_diagnostics;
 
     Diagnostics_sptr step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"test2_full2_per_step.h5"));
     full2_diagnostics.append(step_full2_diag_sptr);
    
    Diagnostics_sptr second_step_full2_diag_sptr(new  Diagnostics_full2(bunch_sptr,"test2_second_full2_per_step.h5"));
    full2_diagnostics.append(second_step_full2_diag_sptr);
     
     Diagnostics_sptr turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"test2_particles_per_turn.h5"));
     particle_diagnostics.append(turn_particles_diag_sptr);
   
     Diagnostics_sptr second_turn_particles_diag_sptr(new  Diagnostics_particles(bunch_sptr,"test2_second_particles_per_turn.h5"));
     particle_diagnostics.append(second_turn_particles_diag_sptr);
    
     
     Bunch_with_diagnostics bunch_with_diagnostics(bunch,full2_diagnostics, particle_diagnostics);
}*/
