#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics.h"
#include "synergia/simulation/diagnostics_normal_form.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/simulation/fast_normal_form.h"



BOOST_GLOBAL_FIXTURE(MPI_fixture);

const int num_macro_particles = 4096;
const double num_real_particles = 1.0e11;
struct Normalform_fixture
{

Normalform_fixture():
    comm_sptr(new Commxx),
    bunch_sptr()
  {
      BOOST_TEST_MESSAGE("Normalform_fixture");
      
      
      lattice_sptr = Lattice_sptr(
          new Lattice(read_lsexpr_file("lattices/nonlinear_lattice.lsx")));
      
      lattice_simulator=Lattice_simulator(lattice_sptr, 3);
       nf_sptr=lattice_simulator.get_normal_form_sptr();
       fnf_sptr=Fast_normal_form_sptr(new Fast_normal_form(*(nf_sptr)));
       bunch_sptr = Bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(),
                                        num_macro_particles, num_real_particles,
                                        comm_sptr));
       
       for (int part = 0; part < bunch_sptr->get_local_num(); ++part) {
                for (int i = 0; i < 6; i += 1) {
                        bunch_sptr->get_local_particles()[part][i] = 10.0 * part + (1.0 + part
                                * part / 1000.0) * i + 1000 * bunch_sptr->get_comm().get_rank();
                }
                bunch_sptr->get_local_particles()[part][Bunch::id] = part + bunch_sptr->get_total_num()
                            * bunch_sptr->get_comm().get_rank();
       }
     
  }

  ~Normalform_fixture()
  {
    BOOST_TEST_MESSAGE("teardown Normalform_fixture");
  }

  Lattice_sptr lattice_sptr;
  Lattice_simulator lattice_simulator;
  Normal_form_sage_sptr nf_sptr;
  Fast_normal_form_sptr  fnf_sptr;  
  Commxx_sptr comm_sptr;
  Bunch_sptr bunch_sptr;

    
};


const double tolerance = 1.0e-11;

BOOST_FIXTURE_TEST_CASE(construct, Normalform_fixture)
{
      Diagnostics_normal_form diag_nf(fnf_sptr,"nf_diagnostics.h5");
    
}

BOOST_FIXTURE_TEST_CASE(write_, Normalform_fixture)
{
     Diagnostics_normal_form diag_nf(fnf_sptr,"nf_diagnostics.h5");
     diag_nf.set_bunch_sptr(bunch_sptr);
     diag_nf.update();
     diag_nf.write();
     
}

// BOOST_FIXTURE_TEST_CASE(write_, Bunch_sptr_fixture)
// {
//     Diagnostics_full2 diagnostics("diagnostics_full2_mpi.h5");
//     diagnostics.set_bunch_sptr(bunch_sptr);
//     diagnostics.update();
//     diagnostics.write();
// }
// 
// BOOST_FIXTURE_TEST_CASE(write_several, Bunch_sptr_fixture)
// {
//     Diagnostics_full2 diagnostics("diagnostics_full2_mpi.h5");
//     diagnostics.set_bunch_sptr(bunch_sptr);
//     diagnostics.update_and_write();
//     diagnostics.update_and_write();
//     diagnostics.update_and_write();
//     diagnostics.update_and_write();
// }
