#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/wake_field.h"
// #include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
// #include "synergia/foundation/physical_constants.h"
// #include "synergia/simulation/step.h"
// #include "synergia/simulation/operator.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture)
const double tolerance = 1.0e-11;
// 
// 
// double const  orbit_length=160.;
// double const bunch_spacing=4.5;
// int const  zgrid=40;
// std::string const pipe_symmetry("x_parallel");
// 
// 
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


BOOST_AUTO_TEST_CASE(test_constructor)
{
Wake_field wakef("test_wake.dat","XLXTYLYTZ");
}

 
BOOST_AUTO_TEST_CASE(wake_reading_3XLYLZ)
 {
 Wake_field wakef("test_wake_cc.dat", "XLYLZ"); // three columns file
 
     
      BOOST_CHECK_CLOSE(wakef.get_z_coord()[0], 0.00614575, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_z_coord()[1], 0.0136406, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_z_coord()[2], 0.0241333, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_z_coord()[3], 0.037624, tolerance);
      
      BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  9.30907e+08, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_z_wake()[1],  5.79371e+08, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_z_wake()[2],  3.78051e+08 , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_z_wake()[3],  2.55723e+08, tolerance);
      
      BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],   -7.40775e+10   , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_xw_lead()[1],   -8.2791e+10, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_xw_lead()[2],   -8.95615e+10 , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_xw_lead()[3],   -9.39278e+10  , tolerance);
      
      BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0],   -7.40775e+10   , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_yw_lead()[1],   -8.2791e+10, tolerance);
      BOOST_CHECK_CLOSE(wakef.get_yw_lead()[2],   -8.95615e+10 , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_yw_lead()[3],   -9.39278e+10  , tolerance);
      
      BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0],   0.   , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_xw_trail()[3],   0., tolerance);
      BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0],   0.   , tolerance);
      BOOST_CHECK_CLOSE(wakef.get_yw_trail()[3],   0., tolerance);
     
 }
 
 BOOST_AUTO_TEST_CASE(wake_reading_4XLYLZ)
 {
 Wake_field wakef("test_wake.dat",  "XLYLZ"); // four columns file
 
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[0], 0.00614575, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[1], 0.0136406, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[2], 0.0241333, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[3], 0.037624, tolerance);
 
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  9.30907e+08, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[1],  5.79371e+08, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[2],  3.78051e+08 , tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[3],  2.55723e+08, tolerance);

BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],   -7.40775e+10   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[1],   -8.2791e+10, tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[3],   -9.39278e+10  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0],   -1.26537e+11   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[1],   -1.4772e+11 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[3],   -1.77089e+11  , tolerance);


BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0],   0.   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_trail()[3],   0., tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0],   0.   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[3],   0., tolerance);
 
 }
 
 BOOST_AUTO_TEST_CASE(wake_reading_4XLXTYLYTZ)
 {
 Wake_field wakef("test_wake.dat",  "XLXTYLYTZ"); // four columns file
 
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[0], 0.00614575, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[1], 0.0136406, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[2], 0.0241333, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[3], 0.037624, tolerance);
 
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  9.30907e+08, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[1],  5.79371e+08, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[2],  3.78051e+08 , tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[3],  2.55723e+08, tolerance);

BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],   -7.40775e+10   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[1],   -8.2791e+10, tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[3],   -9.39278e+10  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0],   -1.26537e+11   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_trail()[1],   -1.4772e+11 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_trail()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_trail()[3],   -1.77089e+11  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0],   -7.40775e+10   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[1],   -8.2791e+10, tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[3],   -9.39278e+10  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0],   -1.26537e+11   , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[1],   -1.4772e+11 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[3],   -1.77089e+11  , tolerance);
 
 }


BOOST_AUTO_TEST_CASE(wake_reading_5XLXTYLYTZ)
{
Wake_field wakef("test_5XLXTYLYT.dat", "XLXTYLYT"); // four columns file
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[0],  -0.534401, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[1],  -0.513238, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[2],  -0.492502, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[3],  -0.472194, tolerance);
 
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  0., tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[1],  0., tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[2],  0. , tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[3],  0., tolerance);
 
 
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],   6.80002e+07   , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_lead()[1],   -8.2791e+10, tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[3],   5.65593e+07  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0],   6.70002e+07   , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_trail()[1],   -1.4772e+11 , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_trail()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_trail()[3],    5.55593e+07  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0],   6.60002e+07  , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_lead()[1],   -8.2791e+10, tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[3],   5.45593e+07  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0],    6.50002e+07   , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_trail()[1],   -1.4772e+11 , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_trail()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[3],    5.35593e+07  , tolerance);


}

BOOST_AUTO_TEST_CASE(wake_reading_6XLXTYLYTZ)
{
Wake_field wakef("test_gwake.dat", "XLXTYLYTZ"); // four columns file
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[0],  -0.534401, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[1],  -0.513238, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[2],  -0.492502, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_coord()[3],  -0.472194, tolerance);
 
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  117362, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[1],  122459, tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[2],  108740 , tolerance);
 BOOST_CHECK_CLOSE(wakef.get_z_wake()[3],  83079.8, tolerance);
 
 
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],   6.80002e+07   , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_lead()[1],   -8.2791e+10, tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_lead()[3],   5.65593e+07  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0],   6.70002e+07   , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_trail()[1],   -1.4772e+11 , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_xw_trail()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_xw_trail()[3],    5.55593e+07  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0],   6.60002e+07  , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_lead()[1],   -8.2791e+10, tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_lead()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_lead()[3],   5.45593e+07  , tolerance);

BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0],    6.50002e+07   , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_trail()[1],   -1.4772e+11 , tolerance);
//BOOST_CHECK_CLOSE(wakef.get_yw_trail()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(wakef.get_yw_trail()[3],    5.35593e+07  , tolerance);


}

BOOST_AUTO_TEST_CASE(wake_reading_pp)
{
Wake_field wakef("test_wake_pp.dat", "XLXTYLYTZpp"); // four columns file
  BOOST_CHECK_CLOSE(wakef.get_z_coord()[0],  0.000106884, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0], -6.79233e+09, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0], 6.79233e+09, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0], -9.95614e+09, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0], -6.79233e+09, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],   3.01394e+08, tolerance);
  int last=wakef.get_z_coord().size()-1;
   BOOST_CHECK_CLOSE(wakef.get_z_coord()[last],   1536.57, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_lead()[last],  -879923, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_trail()[last],  879923 , tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_lead()[last],   -372155 , tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_trail()[last], -879923, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_z_wake()[last],   -40.5817, tolerance);  
 
}

BOOST_AUTO_TEST_CASE(Fwake)
{
Wake_field wakef("Fwake.dat", "XLXTYLYTZpp"); // four columns file
  BOOST_CHECK_EQUAL(wakef.get_z_coord().size(), 3000);

  BOOST_CHECK_CLOSE(wakef.get_z_coord()[0],  -0.534401, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],  1.20964e+08, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0], -1.20964e+08, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0], 8.68572e+08, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0], 1.21505e+08, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  205815, tolerance);
 
  int last=wakef.get_z_coord().size()-1;
  BOOST_CHECK_CLOSE(wakef.get_z_coord()[last],   1859.06, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_lead()[last],   -4.25699e+06, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_trail()[last],  4.25699e+06 , tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_lead()[last],   -3.39949e+06 , tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_trail()[last], -4.25699e+06, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_z_wake()[last],   -34.4432, tolerance);  
 
}
 