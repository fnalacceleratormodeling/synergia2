#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/lattice_simulator.h"
#include "lattice_fixture.h"
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/simulation/fast_normal_form.h"
#include <basic_toolkit/VectorD.h>
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
const double tolerance1 = 1.0e-10;
double stdx=0.012;
double stdy=0.01;
double stdct=3.;


BOOST_FIXTURE_TEST_CASE(stationary_actions_order1,Nonlinearlattice_fixture)
{
    int map_order = 4;
  

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    std::vector<double> ls_sa=lattice_simulator.get_stationary_actions(stdx,stdy,stdct);
   // std::cout<<" lattice sa ="<<ls_sa[0]<<", "<<ls_sa[1]<<", "<<ls_sa[2]<<std::endl;   
     
    Normal_form_sage_sptr nf_sptr(lattice_simulator.get_normal_form_sptr());   
    std::vector<double> chef_sa=nf_sptr->stationaryActions(stdx,stdy,stdct/pconstants::c);
   // std::cout<<" chef sa ="<<chef_sa[0]<<", "<<chef_sa[1]<<", "<<chef_sa[2]<<std::endl;   

     
   Fast_normal_form  fnf(*(nf_sptr));
   fnf.get_stationary_actions(stdx,stdy,stdct);
   std::vector<double> fast_sa=fnf.get_stationary_actions(stdx,stdy,stdct);
  // std::cout<<" fast sa ="<<fast_sa[0]<<", "<<fast_sa[1]<<", "<<fast_sa[2]<<std::endl;   

       
    BOOST_CHECK_CLOSE(chef_sa[0] ,  fast_sa[0] ,  tolerance);
    BOOST_CHECK_CLOSE(chef_sa[1] ,  fast_sa[1] ,  tolerance);
    BOOST_CHECK_CLOSE(chef_sa[2] ,  fast_sa[2] ,  tolerance);
    BOOST_CHECK_CLOSE(ls_sa[0] ,  fast_sa[0] ,  tolerance);
    BOOST_CHECK_CLOSE(ls_sa[1] ,  fast_sa[1] ,  tolerance);
    BOOST_CHECK_CLOSE(ls_sa[2] ,  fast_sa[2] ,  tolerance);
    
}



BOOST_FIXTURE_TEST_CASE(convert_to_normal_and_back, Nonlinearlattice_fixture)

{
    int map_order = 4;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

   // MArray1d clo=lattice_simulator.get_closed_orbit();
     //std::cout<<" closed orbit ="<<clo[0]<<", "<<clo[1]<<", "<<clo[2]<<", "<<clo[3]<<", "<<clo[4]<<", "<<clo[5]<<std::endl;


    Normal_form_sage_sptr nf_sptr(lattice_simulator.get_normal_form_sptr());   
    Fast_normal_form  fnf(*(nf_sptr));


   // const double test_points[] = { 2.0e-3, 7.76e-6, 1.0e-3, 1.86e-5, 6.0e-1, 1.0e-4 };
    const double test_points[] = { 2.0e-3, 7.76e-6, 1.0e-3, 1.86e-5, 6.0e-3, 1.0e-5 };

    
         
    MArray1d xyz_coord(boost::extents[6]);
    xyz_coord[0]=test_points[0];
    xyz_coord[1]=test_points[1];
    xyz_coord[2]=test_points[2];
    xyz_coord[3]=test_points[3];
    xyz_coord[4]=test_points[4];
    xyz_coord[5]=test_points[5];
    MArray1dc nf_coord(boost::extents[6]);      
    fnf.convert_to_normal_form(xyz_coord,nf_coord);
    
//      std::cout<<"  nf_fast="<< nf_coord[0]<<",  "
//      << nf_coord[1]<<",  "
//      << nf_coord[2]<<",  "
//      << nf_coord[3]<<",  "
//      << nf_coord[4]<<",  "
//      << nf_coord[5]<<",  "
//      <<std::endl;
//     
    
     Vector xyz_chef(6);
     for (int i=0;i<6;++i){
       int chef_index=get_chef_index(i);
       xyz_chef(chef_index)=xyz_coord[i];
     }
     VectorC nf_chef(6);
 /*    
       std::cout<<"  xyz_chef start="
      << xyz_chef(0)<<",  "
      << xyz_chef(1)<<",  "
      << xyz_chef(2)<<",  "
      << xyz_chef(3)<<",  "
      << xyz_chef(4)<<",  "
      << xyz_chef(5)
      <<std::endl;
     */
     
     nf_sptr-> cnvDataToNormalForm(xyz_chef,  nf_chef);

 /*    
      std::cout<<"  nf_chef="
//      << nf_chef(get_chef_index(0))<<",  "
//      << nf_chef(get_chef_index(1))<<",  "
//      << nf_chef(get_chef_index(2))<<",  "
//      << nf_chef(get_chef_index(3))<<",  "
//      << nf_chef(get_chef_index(4))<<",  "
//      << nf_chef(get_chef_index(5))<<",  "
     << nf_chef(0)<<",  "
     << nf_chef(1)<<",  "
     << nf_chef(2)<<",  "
     << nf_chef(3)<<",  "
     << nf_chef(4)<<",  "
     << nf_chef(5)
     <<std::endl;*/

     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(0)].real() ,  nf_coord[0].real() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(1)].real() ,  nf_coord[1].real() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(2)].real() ,  nf_coord[2].real() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(3)].real() ,  nf_coord[3].real() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(4)].real() ,  nf_coord[4].real() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(5)].real() ,  nf_coord[5].real() ,  tolerance1);
     
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(0)].imag() ,  nf_coord[0].imag() ,  tolerance1);      
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(1)].imag() ,  nf_coord[1].imag() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(2)].imag() ,  nf_coord[2].imag() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(3)].imag() ,  nf_coord[3].imag() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(4)].imag() ,  nf_coord[4].imag() ,  tolerance1);
     BOOST_CHECK_CLOSE(nf_chef[get_chef_index(5)].imag() ,  nf_coord[5].imag() ,  tolerance1);
     
     
     
     MArray1d xyz_back(boost::extents[6]);
     fnf.convert_from_normal_form(nf_coord, xyz_back);
/*
     std::cout<<"  fast_xyz_back="
     << xyz_back[0]<<",  "
     << xyz_back[1]<<",  "
     << xyz_back[2]<<",  "
     << xyz_back[3]<<",  "
     << xyz_back[4]<<",  "
     << xyz_back[5]
     <<std::endl;
    */
     
     
     nf_sptr->cnvDataFromNormalForm(nf_chef,  xyz_chef);
//           std::cout<<"  xyz_chef_back="
//             
// //      << xyz_chef(get_chef_index(0))<<",  "
// //      << xyz_chef(get_chef_index(1))<<",  "
// //      << xyz_chef(get_chef_index(2))<<",  "
// //      << xyz_chef(get_chef_index(3))<<",  "
// //      << xyz_chef(get_chef_index(4))<<",  "
// //      << xyz_chef(get_chef_index(5))  
//      << xyz_chef(0)<<",  "
//      << xyz_chef(1)<<",  "
//      << xyz_chef(2)<<",  "
//      << xyz_chef(3)<<",  "
//      << xyz_chef(4)<<",  "
//      << xyz_chef(5)
//      <<std::endl;
//      
     
     BOOST_CHECK_CLOSE(xyz_chef[get_chef_index(0)] ,  xyz_back[0] ,  tolerance1);
     BOOST_CHECK_CLOSE(xyz_chef[get_chef_index(1)] ,  xyz_back[1] ,  tolerance1);
     BOOST_CHECK_CLOSE(xyz_chef[get_chef_index(2)] ,  xyz_back[2] ,  tolerance1);
     BOOST_CHECK_CLOSE(xyz_chef[get_chef_index(3)] ,  xyz_back[3] ,  tolerance1);
     BOOST_CHECK_CLOSE(xyz_chef[get_chef_index(4)] ,  xyz_back[4] ,  tolerance1);
     BOOST_CHECK_CLOSE(xyz_chef[get_chef_index(5)] ,  xyz_back[5] ,  tolerance1);
     
     
     
    
}




BOOST_FIXTURE_TEST_CASE(serialize_, Nonlinearlattice_fixture)

{
    int map_order = 4;
    const double test_points[] = { 2.0e-3, 7.76e-6, 1.0e-3, 1.86e-5, 6.0e-3, 1.0e-5 };
    MArray1d xyz_coord(boost::extents[6]);
    xyz_coord[0]=test_points[0];
    xyz_coord[1]=test_points[1];
    xyz_coord[2]=test_points[2];
    xyz_coord[3]=test_points[3];
    xyz_coord[4]=test_points[4];
    xyz_coord[5]=test_points[5];
    
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

   // MArray1d clo=lattice_simulator.get_closed_orbit();
     //std::cout<<" closed orbit ="<<clo[0]<<", "<<clo[1]<<", "<<clo[2]<<", "<<clo[3]<<", "<<clo[4]<<", "<<clo[5]<<std::endl;


    Normal_form_sage_sptr nf_sptr(lattice_simulator.get_normal_form_sptr());   
    Fast_normal_form  fnf(*(nf_sptr));      
    std::vector<double> fast_sa=fnf.get_stationary_actions(stdx,stdy,stdct);
    
    MArray1dc nf_coord(boost::extents[6]);      
    fnf.convert_to_normal_form(xyz_coord,nf_coord);
    
     MArray1d xyz_back(boost::extents[6]);
     fnf.convert_from_normal_form(nf_coord, xyz_back);
    
    
    xml_save(fnf, "fast_normal_form.xml");                
    Fast_normal_form   fnf_load;
    xml_load(fnf_load, "fast_normal_form.xml");
    
    
    
    std::vector<double> load_sa=fnf_load.get_stationary_actions(stdx,stdy,stdct);
    
    BOOST_CHECK_CLOSE(load_sa[0] ,  fast_sa[0] ,  tolerance);
    BOOST_CHECK_CLOSE(load_sa[1] ,  fast_sa[1] ,  tolerance);
    BOOST_CHECK_CLOSE(load_sa[2] ,  fast_sa[2] ,  tolerance);
    
    MArray1dc nf_loaded(boost::extents[6]); 
    fnf_load.convert_to_normal_form(xyz_coord,nf_loaded);
  
    BOOST_CHECK_CLOSE(nf_loaded[0].real() ,  nf_coord[0].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[1].real() ,  nf_coord[1].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[2].real() ,  nf_coord[2].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[3].real() ,  nf_coord[3].real() ,  tolerance);     
    BOOST_CHECK_CLOSE(nf_loaded[4].real() ,  nf_coord[4].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[5].real() ,  nf_coord[5].real() ,  tolerance);

    BOOST_CHECK_CLOSE(nf_loaded[0].imag() ,  nf_coord[0].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[1].imag() ,  nf_coord[1].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[2].imag() ,  nf_coord[2].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[3].imag() ,  nf_coord[3].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[4].imag() ,  nf_coord[4].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_loaded[5].imag() ,  nf_coord[5].imag() ,  tolerance);
    
    MArray1d loaded_back(boost::extents[6]);
    fnf_load.convert_from_normal_form(nf_loaded,loaded_back );
    
    BOOST_CHECK_CLOSE(loaded_back[0] ,  xyz_back[0] ,  tolerance);
    BOOST_CHECK_CLOSE(loaded_back[1] ,  xyz_back[1] ,  tolerance);
    BOOST_CHECK_CLOSE(loaded_back[2] ,  xyz_back[2] ,  tolerance);
    BOOST_CHECK_CLOSE(loaded_back[3] ,  xyz_back[3] ,  tolerance);
    BOOST_CHECK_CLOSE(loaded_back[4] ,  xyz_back[4] ,  tolerance);
    BOOST_CHECK_CLOSE(loaded_back[5] ,  xyz_back[5] ,  tolerance);
    
    
    
    archive_save<Fast_normal_form, boost::archive::text_oarchive > (fnf, "fast_normal_form.txt",false);     
    Fast_normal_form   fnf_text;
    archive_load<Fast_normal_form, boost::archive::text_iarchive > (fnf_text, "fast_normal_form.txt"); 
    
    
    
    std::vector<double> text_sa=fnf_text.get_stationary_actions(stdx,stdy,stdct);
    BOOST_CHECK_CLOSE(text_sa[0] ,  fast_sa[0] ,  tolerance);
    BOOST_CHECK_CLOSE(text_sa[1] ,  fast_sa[1] ,  tolerance);
    BOOST_CHECK_CLOSE(text_sa[2] ,  fast_sa[2] ,  tolerance);
    
    MArray1dc nf_text(boost::extents[6]); 
    fnf_text.convert_to_normal_form(xyz_coord,nf_text);
    BOOST_CHECK_CLOSE(nf_text[0].real() ,  nf_coord[0].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[1].real() ,  nf_coord[1].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[2].real() ,  nf_coord[2].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[3].real() ,  nf_coord[3].real() ,  tolerance);     
    BOOST_CHECK_CLOSE(nf_text[4].real() ,  nf_coord[4].real() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[5].real() ,  nf_coord[5].real() ,  tolerance);

    BOOST_CHECK_CLOSE(nf_text[0].imag() ,  nf_coord[0].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[1].imag() ,  nf_coord[1].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[2].imag() ,  nf_coord[2].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[3].imag() ,  nf_coord[3].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[4].imag() ,  nf_coord[4].imag() ,  tolerance);
    BOOST_CHECK_CLOSE(nf_text[5].imag() ,  nf_coord[5].imag() ,  tolerance);
    
   
     
    MArray1d text_back(boost::extents[6]);
    fnf_text.convert_from_normal_form(nf_text,text_back );
    BOOST_CHECK_CLOSE(text_back[0] ,  xyz_back[0] ,  tolerance);
    BOOST_CHECK_CLOSE(text_back[1] ,  xyz_back[1] ,  tolerance);
    BOOST_CHECK_CLOSE(text_back[2] ,  xyz_back[2] ,  tolerance);
    BOOST_CHECK_CLOSE(text_back[3] ,  xyz_back[3] ,  tolerance);
    BOOST_CHECK_CLOSE(text_back[4] ,  xyz_back[4] ,  tolerance);
    BOOST_CHECK_CLOSE(text_back[5] ,  xyz_back[5] ,  tolerance);   
}



BOOST_FIXTURE_TEST_CASE(convert_coordinates, Nonlinearlattice_fixture)

{
    int map_order = 4;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);



    Normal_form_sage_sptr nf_sptr(lattice_simulator.get_normal_form_sptr());   
    Fast_normal_form  fnf(*(nf_sptr));

    const double test_points[] = { 2.0e-3, 7.76e-6, 1.0e-3, 1.86e-5, 6.0e-3, 1.0e-5 };

    
    int npoints=8;     
    MArray2d coords(boost::extents[npoints][7]);     
    for (int i=0;i<npoints;++i){
        double fact=-0.3+i*0.1;
        coords[i][0]=fact*test_points[0];
        coords[i][1]=fact*test_points[1];
        coords[i][2]=fact*test_points[2];
        coords[i][3]=fact*test_points[3];
        coords[i][4]=fact*test_points[4];
        coords[i][5]=fact*test_points[5];
        coords[i][6]=i;
    }  
    MArray2d fast_coords(coords); 

    
     
    lattice_simulator.convert_xyz_to_normal(coords);
    fnf.convert_xyz_to_normal(fast_coords);

    for (int i=0;i<npoints;++i){
        for (int j=0;j<7;++j){
            BOOST_CHECK_CLOSE( coords[i][j] , fast_coords[i][j],  tolerance);
        }
//          std::cout<<" coords["<<i<<"][0]="<< coords[i][0]<< std::endl;
//          std::cout<<" coords["<<i<<"][1]="<< coords[i][1]<< std::endl;
//          std::cout<<" coords["<<i<<"][2]="<< coords[i][2]<< std::endl;
//          std::cout<<" coords["<<i<<"][3]="<< coords[i][3]<< std::endl;
//          std::cout<<" coords["<<i<<"][4]="<< coords[i][4]<< std::endl;
//          std::cout<<" coords["<<i<<"][5]="<< coords[i][5]<< std::endl;
//          std::cout<<" coords["<<i<<"][6]="<< coords[i][6]<< std::endl;
//          std::cout<< std::endl;
//          
//          
//          std::cout<<" fcoords["<<i<<"][0]="<< fast_coords[i][0]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][1]="<< fast_coords[i][1]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][2]="<< fast_coords[i][2]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][3]="<< fast_coords[i][3]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][4]="<< fast_coords[i][4]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][5]="<< fast_coords[i][5]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][6]="<< fast_coords[i][6]<< std::endl;
//          std::cout<< std::endl;
//          std::cout<< std::endl;
//                   
     }
     
    lattice_simulator.convert_normal_to_xyz(coords); 
    fnf.convert_normal_to_xyz(fast_coords);
    for (int i=0;i<npoints;++i){
         for (int j=0;j<7;++j){
            BOOST_CHECK_CLOSE( coords[i][j] , fast_coords[i][j],  tolerance1);
        }
        
//          std::cout<<" coords["<<i<<"][0]="<< coords[i][0]<< std::endl;
//          std::cout<<" coords["<<i<<"][1]="<< coords[i][1]<< std::endl;
//          std::cout<<" coords["<<i<<"][2]="<< coords[i][2]<< std::endl;
//          std::cout<<" coords["<<i<<"][3]="<< coords[i][3]<< std::endl;
//          std::cout<<" coords["<<i<<"][4]="<< coords[i][4]<< std::endl;
//          std::cout<<" coords["<<i<<"][5]="<< coords[i][5]<< std::endl;
//          std::cout<<" coords["<<i<<"][6]="<< coords[i][6]<< std::endl;
//          std::cout<< std::endl; 
//                   
//          std::cout<<" fcoords["<<i<<"][0]="<< fast_coords[i][0]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][1]="<< fast_coords[i][1]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][2]="<< fast_coords[i][2]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][3]="<< fast_coords[i][3]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][4]="<< fast_coords[i][4]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][5]="<< fast_coords[i][5]<< std::endl;
//          std::cout<<" fcoords["<<i<<"][6]="<< fast_coords[i][6]<< std::endl;
//          std::cout<< std::endl;
//          std::cout<< std::endl;
                   
        
    }
    
}  
     
  