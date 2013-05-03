#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/wake_field.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/multi_array_check_equal.h"


BOOST_GLOBAL_FIXTURE(MPI_fixture)
const double tolerance = 1.0e-11;


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

BOOST_AUTO_TEST_CASE(Fwake_multiply)
{
  Wake_field wakef("Fwake.dat", "XLXTYLYTZpp"); // four columns file
  BOOST_CHECK_EQUAL(wakef.get_z_coord().size(), 3000);
 
  double mltp=2.4;
  wakef.multiply_xw_lead(mltp);
  wakef.multiply_xw_trail(mltp);
  wakef.multiply_yw_lead(mltp);
  wakef.multiply_yw_trail(mltp);
  wakef.multiply_z_wake(mltp);

  BOOST_CHECK_CLOSE(wakef.get_z_coord()[0],  -0.534401, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_lead()[0],  1.20964e+08*mltp, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_trail()[0], -1.20964e+08*mltp, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_lead()[0], 8.68572e+08*mltp, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_trail()[0], 1.21505e+08*mltp, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_z_wake()[0],  205815*mltp, tolerance);
 
  int last=wakef.get_z_coord().size()-1;
  BOOST_CHECK_CLOSE(wakef.get_z_coord()[last],   1859.06, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_lead()[last],   -4.25699e+06*mltp, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_xw_trail()[last],  4.25699e+06*mltp , tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_lead()[last],   -3.39949e+06*mltp , tolerance);
  BOOST_CHECK_CLOSE(wakef.get_yw_trail()[last], -4.25699e+06*mltp, tolerance);
  BOOST_CHECK_CLOSE(wakef.get_z_wake()[last],   -34.4432*mltp, tolerance);  
 
}


BOOST_AUTO_TEST_CASE(serialize_)
{
  Wake_field wakef("Fwake.dat", "XLXTYLYTZpp"); // four columns file
  xml_save(wakef, "wakef.xml");
 
  Wake_field wake_load;
  xml_load(wake_load, "wakef.xml");
   
   
  
   BOOST_CHECK_EQUAL(wakef.get_wake_type(), wake_load.get_wake_type());
   BOOST_CHECK_EQUAL(wakef.get_wake_file_name(), wake_load.get_wake_file_name());
   BOOST_CHECK_EQUAL(wakef.get_istart() , wake_load.get_istart() );
   BOOST_CHECK_EQUAL(wakef.get_zstart() , wake_load.get_zstart() );
   BOOST_CHECK_EQUAL(wakef.get_delta_z() , wake_load.get_delta_z() );
   multi_array_check_equal(wakef.get_z_coord(), wake_load.get_z_coord(),  tolerance);
   multi_array_check_equal(wakef.get_xw_lead() , wake_load.get_xw_lead() ,  tolerance);
   multi_array_check_equal(wakef.get_xw_trail() , wake_load.get_xw_trail() ,  tolerance);
   multi_array_check_equal(wakef.get_yw_lead() , wake_load.get_yw_lead() ,  tolerance); 
   multi_array_check_equal(wakef.get_yw_trail() , wake_load.get_yw_trail() ,  tolerance); 
   multi_array_check_equal(wakef.get_z_wake() , wake_load.get_z_wake() ,  tolerance); 
        
   
}

 
 
  



 
