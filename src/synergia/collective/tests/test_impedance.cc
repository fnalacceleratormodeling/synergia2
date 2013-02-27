#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/impedance.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/operator.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

 const double tolerance = 1.0e-11;


double const  orbit_length=160.;
double const bunch_spacing=5.;
int const  zgrid=40;



// const double mass = 100.0;
// const double total_energy = 125.0;
// const int total_num = 100;
// const double real_num = 2.0e12;
// const double step_length = 1.23;
// struct Fixture
// {
//     Fixture() :
//         four_momentum(mass, total_energy), reference_particle(
//                 pconstants::proton_charge, four_momentum),
//                 comm_sptr(new Commxx),
//                  bunch(reference_particle, total_num,
//                         real_num, comm_sptr), step(step_length)
//     {
//         BOOST_TEST_MESSAGE("setup fixture");
//     }
//     ~Fixture()
//     {
//         BOOST_TEST_MESSAGE("teardown fixture");
//     }
// 
//     Four_momentum four_momentum;
//     Reference_particle reference_particle;
//     Commxx_sptr comm_sptr;
//     Bunch bunch;
//     Step step;
// };
// 
// void
// dummy_populate(Bunch &bunch, int offset = 0)
// {
//     for (int part = 0; part < bunch.get_local_num(); ++part) {
//         // coordinates
//         for (int i = 0; i < 6; i += 2) {
//             bunch.get_local_particles()[part][i] = 10.0 * (part + offset) + i;
//         }
//         // momenta
//         for (int i = 1; i < 6; i += 2) {
//             bunch.get_local_particles()[part][i] = 1e-4 * (10.0 * (part
//                     + offset) + i);
//         }
//     }
// }
// 
// 
// 
// 
 BOOST_AUTO_TEST_CASE(test_constructor)
 {
 Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid, orbit_length, bunch_spacing,10);
 }
 
BOOST_AUTO_TEST_CASE(test_copy)
 {
 Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid, orbit_length, bunch_spacing,10);
 Impedance imped1(imped);
 imped1.set_z_grid(30);
 
//  std::cout<<" imped zgrid="<<imped.get_z_grid()<<std::endl;
//  std::cout<<" imped wake file="<<imped.get_wake_field_sptr()->get_wake_file_name()<<std::endl;
//  
//  std::cout<<" imped1 zgrid="<<imped1.get_z_grid()<<std::endl;
//  std::cout<<" imped1 wake file="<<imped1.get_wake_field_sptr()->get_wake_file_name()<<std::endl;
 
 Impedance_sptr imped_sptr=Impedance_sptr(new Impedance("test_wake_pp.dat", "XLXTYLYTZpp",
							zgrid, orbit_length, bunch_spacing,100));
 Impedance imped2(*imped_sptr);
 } 


// 
// BOOST_FIXTURE_TEST_CASE(test_apply, Fixture)
//   {
//     Bunch bunch(reference_particle, total_num, real_num, comm_sptr);
//     dummy_populate(bunch);
//     Step step(step_length);
//     double time_step=10.;
//     Impedance imped("test_wake.dat", orbit_length, bunch_spacing, zgrid, "circular",10); // four columns file
//     const int verbosity = 4;
//     Logger logger(0);
//     imped.apply(bunch, time_step, step, verbosity, logger);
// }
