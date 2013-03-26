#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/hdf5_file.h"
#include "gaussian_charge_density.h"
#include "space_charge_bunch_fixtures.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/commxx_divider.h"

 BOOST_GLOBAL_FIXTURE(MPI_fixture)

 const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 32;
    grid_shape[1] = 16;
    grid_shape[2] = 18;
    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.2;
    size[2]=0.05;
    Commxx_sptr comm_f_sptr(new Commxx);
    Space_charge_rectangular space_charge(comm_f_sptr,size, grid_shape);

}

BOOST_FIXTURE_TEST_CASE(get_local_charge_density, Ellipsoidal_bunch_fixture)
{

    std::vector<int > grid_shape(3);
    grid_shape[0] = 32;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.2;
    size[2]=0.05;
    Space_charge_rectangular space_charge(size, grid_shape);
    space_charge.set_fftw_helper(bunch.get_comm_sptr(), false);
    Rectangular_grid_sptr rho(space_charge.get_charge_density(bunch));

  //double local_charge=bunch.get_local_num()*real_num*pconstants::e/total_num;
  // BOOST_CHECK_CLOSE(local_rho->get_normalization(), local_charge, 100*tolerance);

}

BOOST_FIXTURE_TEST_CASE(get_charge_density1, Toy_bunch_fixture_xyz)
{

    Space_charge_rectangular space_charge(physical_size, grid_shape);

    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;
    Rectangular_grid_sptr rho(space_charge.get_charge_density(bunch));
    std::vector<int > center(3);
    center[0] = grid_shape[0] / 2;
    center[1] = grid_shape[1] / 2;
    center[2] = grid_shape[2] / 2;

    int size_com=bunch.get_comm().get_size();

    for (int i = center[0] - 1; i < center[0] + 1; ++i) {
        for (int j = center[1] - 1; j < center[1] + 1; ++j) {
            for (int k = center[2] - 1; k < center[2] + 1; ++k) {
                expected[i][j][k] = 0.125 * density_norm*size_com;
            }
        }
    }
    multi_array_check_equal(rho->get_grid_points(), expected,
            100 * tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_phi_local, Ellipsoidal_bunch_fixture)
{
    //grid_shape[0] = 16;
    //grid_shape[1] = 24;
    //grid_shape[2] = 32;

    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.2;
    size[2]=0.05;

    Space_charge_rectangular space_charge(size, grid_shape);
    Rectangular_grid_sptr rho(space_charge.get_charge_density(bunch));

    space_charge.set_fftw_helper(bunch.get_comm_sptr(), false);
    Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho));
    
 }

BOOST_FIXTURE_TEST_CASE(get_phi_local1, Ellipsoidal_bunch_fixture)
{
    //grid_shape[0] = 16;
    //grid_shape[1] = 24;
    //grid_shape[2] = 32;

    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.2;
    size[2]=0.05;

    Space_charge_rectangular space_charge(size, grid_shape);
    Rectangular_grid_sptr rho(space_charge.get_charge_density(bunch));
    
    int optimal_number=3;
    Commxx_sptr comm_spc=make_optimal_spc_comm(bunch.get_comm_sptr(), optimal_number);

    space_charge.set_fftw_helper(comm_spc, false);
    if (space_charge.get_comm_sptr()->has_this_rank())
	Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho));
    
 }
 
 
 BOOST_FIXTURE_TEST_CASE(get_phi_local2, Ellipsoidal_bunch_fixture)
{
    //grid_shape[0] = 16;
    //grid_shape[1] = 24;
    //grid_shape[2] = 32;

    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.2;
    size[2]=0.05;

    Space_charge_rectangular space_charge(size, grid_shape);
    Rectangular_grid_sptr rho(space_charge.get_charge_density(bunch));
    
    int optimal_number=2; 
    if (bunch.get_comm_sptr()->get_size() %  optimal_number==0) {
      Commxx_divider comm_divider(optimal_number, false);
      Commxx_sptr comm_spc=comm_divider.get_commxx_sptr(bunch.get_comm_sptr());
      space_charge.set_fftw_helper(comm_spc, true);
      Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho));
    }
    
 }

   BOOST_FIXTURE_TEST_CASE(get_phi_local3, Ellipsoidal_bunch_fixture)
  {

     std::vector<int > grid_shape1(3);
     grid_shape1[0] =32;
     grid_shape1[1] = 16;
     grid_shape1[2] = 12;
     std::vector<double > size(3);
     size[0]=0.1;
     size[1]=0.1;
     size[2]=0.1;

     Space_charge_rectangular space_charge(size, grid_shape1);
     std::vector<double > hi(space_charge.get_domain().get_cell_size());
     std::vector<int > shape(space_charge.get_domain().get_grid_shape());

     int nt=18;
     int mt=1;
     int pt=2;


    Rectangular_grid_sptr rho_grid(new Rectangular_grid(space_charge.get_domain_sptr()));

    MArray3d_ref rho(rho_grid->get_grid_points());
    for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                for (int k = 0; k < shape[2]; ++k) {
                   rho[i][j][k] =  sin(nt*(i+0.5)*pi/double(shape[0]))*sin(mt*(j+0.5)*pi/double(shape[1]))*cos((2.*pi*pt*k)/double(shape[2]));

                }
            }
        }


//     int size_com=bunch.get_comm().get_size();
//     int lrank=bunch.get_comm().get_rank();
//     space_charge.set_fftw_helper(bunch.get_comm_sptr());
    
     int optimal_number=3;
     Commxx_sptr comm_spc=make_optimal_spc_comm(bunch.get_comm_sptr(), optimal_number);
     space_charge.set_fftw_helper(comm_spc,false);
    
   
    
    if (space_charge.get_comm_sptr()->has_this_rank()){
	int size_com=comm_spc->get_size();    
	int lrank=comm_spc->get_rank();
	Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho_grid));
	std::vector<int> local_shape(phi_local->get_domain().get_grid_shape());

	// std::cout<<" local_shape="<<local_shape[0]<<" rank="<<lrank<<std::endl;

	double xnt=pi*nt/size[0];
	double ymt=pi*mt/size[1];
	double zpt=2.*pi*pt/size[2];
	double Cnmp=-(xnt*xnt+ymt*ymt+zpt*zpt);


    //    local check
	int offset=phi_local->get_lower();
	MArray3d  rhocheck_local(boost::extents[local_shape[0]][shape[1]][shape[2]]);
	MArray3d  phicheck_local(boost::extents[local_shape[0]][shape[1]][shape[2]]);
	for (int i = 0; i < local_shape[0]; ++i) {
	      int i_local=i+offset;
		for (int j = 0; j < shape[1]; ++j) {
		    for (int k = 0; k < shape[2]; ++k) {
		      rhocheck_local[i][j][k] = rho[i_local][j][k];
		      phicheck_local[i][j][k] = -Cnmp* phi_local->get_grid_points()[i][j][k]*phi_local->get_normalization()*epsilon0;
		    }
		}
	    }

	  multi_array_check_equal(rhocheck_local, phicheck_local, 10.e-8);

	std::vector<int> uppers(phi_local->get_uppers());
	std::vector<int> receive_counts(phi_local->get_lengths()), receive_offsets(size_com);
	for (int i=0; i< size_com; ++i) {
	    receive_offsets.at(i) = uppers.at(i)*shape[1]*shape[2]-receive_counts.at(i);
	}

	Rectangular_grid_sptr phi_global_grid(new Rectangular_grid(space_charge.get_domain_sptr()));
	MArray3d_ref phi_global(phi_global_grid->get_grid_points());



	int error =MPI_Allgatherv(reinterpret_cast<void*>(phi_local->get_grid_points().origin()),
		  receive_counts[lrank], MPI_DOUBLE,
		  reinterpret_cast<void*>(phi_global.origin()),
					  &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, 
				                 space_charge.get_comm_sptr()->get());
	if (error != MPI_SUCCESS) {
	    throw std::runtime_error(
		    "MPI error in test_space_charge_rectangular(MPI_Allgatherv in test get_phi_local1)");
	}

	MArray3d  phicheck(boost::extents[shape[0]][shape[1]][shape[2]]);
	for (int i = 0; i < shape[0]; ++i) {
		for (int j = 0; j < shape[1]; ++j) {
		    for (int k = 0; k < shape[2]; ++k) {
			phicheck[i][j][k] = -Cnmp* phi_global[i][j][k]*phi_local->get_normalization()*epsilon0;
		    }
		}
	    }

	multi_array_check_equal(rho, phicheck, 10.e-8);
    }//comm_sptr()->has_this_rank
 }

  BOOST_FIXTURE_TEST_CASE(get_phi_local4, Ellipsoidal_bunch_fixture)
  {

     std::vector<int > grid_shape1(3);
     grid_shape1[0] =32;
     grid_shape1[1] = 16;
     grid_shape1[2] = 12;
     std::vector<double > size(3);
     size[0]=0.1;
     size[1]=0.1;
     size[2]=0.1;

     Space_charge_rectangular space_charge(size, grid_shape1);
     std::vector<double > hi(space_charge.get_domain().get_cell_size());
     std::vector<int > shape(space_charge.get_domain().get_grid_shape());

     int nt=18;
     int mt=1;
     int pt=2;


    Rectangular_grid_sptr rho_grid(new Rectangular_grid(space_charge.get_domain_sptr()));

    MArray3d_ref rho(rho_grid->get_grid_points());
    for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                for (int k = 0; k < shape[2]; ++k) {
                   rho[i][j][k] =  sin(nt*(i+0.5)*pi/double(shape[0]))*sin(mt*(j+0.5)*pi/double(shape[1]))*cos((2.*pi*pt*k)/double(shape[2]));

                }
            }
        }


//     int size_com=bunch.get_comm().get_size();
//     int lrank=bunch.get_comm().get_rank();
//     space_charge.set_fftw_helper(bunch.get_comm_sptr());
    int optimal_number=2;
    if (bunch.get_comm_sptr()->get_size() %  optimal_number==0) {
	Commxx_divider comm_divider(optimal_number, false);
	Commxx_sptr comm_spc=comm_divider.get_commxx_sptr(bunch.get_comm_sptr());
	space_charge.set_fftw_helper(comm_spc, true);
     
    
   
    
   
	int size_com=comm_spc->get_size();    
	int lrank=comm_spc->get_rank();
	Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho_grid));
	std::vector<int> local_shape(phi_local->get_domain().get_grid_shape());

	// std::cout<<" local_shape="<<local_shape[0]<<" rank="<<lrank<<std::endl;

	double xnt=pi*nt/size[0];
	double ymt=pi*mt/size[1];
	double zpt=2.*pi*pt/size[2];
	double Cnmp=-(xnt*xnt+ymt*ymt+zpt*zpt);


    //    local check
	int offset=phi_local->get_lower();
	MArray3d  rhocheck_local(boost::extents[local_shape[0]][shape[1]][shape[2]]);
	MArray3d  phicheck_local(boost::extents[local_shape[0]][shape[1]][shape[2]]);
	for (int i = 0; i < local_shape[0]; ++i) {
	      int i_local=i+offset;
		for (int j = 0; j < shape[1]; ++j) {
		    for (int k = 0; k < shape[2]; ++k) {
		      rhocheck_local[i][j][k] = rho[i_local][j][k];
		      phicheck_local[i][j][k] = -Cnmp* phi_local->get_grid_points()[i][j][k]*phi_local->get_normalization()*epsilon0;
		    }
		}
	    }

	  multi_array_check_equal(rhocheck_local, phicheck_local, 10.e-8);

	std::vector<int> uppers(phi_local->get_uppers());
	std::vector<int> receive_counts(phi_local->get_lengths()), receive_offsets(size_com);
	for (int i=0; i< size_com; ++i) {
	    receive_offsets.at(i) = uppers.at(i)*shape[1]*shape[2]-receive_counts.at(i);
	}

	Rectangular_grid_sptr phi_global_grid(new Rectangular_grid(space_charge.get_domain_sptr()));
	MArray3d_ref phi_global(phi_global_grid->get_grid_points());



	int error =MPI_Allgatherv(reinterpret_cast<void*>(phi_local->get_grid_points().origin()),
		  receive_counts[lrank], MPI_DOUBLE,
		  reinterpret_cast<void*>(phi_global.origin()),
					  &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, 
				                 space_charge.get_comm_sptr()->get());
	if (error != MPI_SUCCESS) {
	    throw std::runtime_error(
		    "MPI error in test_space_charge_rectangular(MPI_Allgatherv in test get_phi_local1)");
	}
    

	MArray3d  phicheck(boost::extents[shape[0]][shape[1]][shape[2]]);
	for (int i = 0; i < shape[0]; ++i) {
		for (int j = 0; j < shape[1]; ++j) {
		    for (int k = 0; k < shape[2]; ++k) {
			phicheck[i][j][k] = -Cnmp* phi_global[i][j][k]*phi_local->get_normalization()*epsilon0;
		    }
		}
	    }

	multi_array_check_equal(rho, phicheck, 10.e-8);
    }//bunch.get_comm_sptr()->get_size() %  optimal_number==0	
   
 }





 BOOST_FIXTURE_TEST_CASE(get_En, Ellipsoidal_bunch_fixture)
{

    std::vector<int > grid_shape1(3);
    grid_shape1[0] =32;
    grid_shape1[1] = 12;
    grid_shape1[2] = 10;
    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.1;
    size[2]=0.1;

    Space_charge_rectangular space_charge(size, grid_shape1);

    int nt=6;
    int mt=3;
    int pt=4;


    Rectangular_grid_sptr rho_grid(new Rectangular_grid(space_charge.get_domain_sptr()));
    std::vector<int > shape(space_charge.get_domain().get_grid_shape());

    MArray3d_ref rho(rho_grid->get_grid_points());
    for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                for (int k = 0; k < shape[2]; ++k) {
                   rho[i][j][k] = sin(nt*(i+0.5)*pi/double(shape[0]))*sin(mt*(j+0.5)*pi/double(shape[1])) *cos((2.*pi*pt*k)/double(shape[2]));
                }
            }
        }

   

   std::vector<double > hi(space_charge.get_domain().get_cell_size());
   int optimal_number=3;
   Commxx_sptr comm_spc=make_optimal_spc_comm(bunch.get_comm_sptr(), optimal_number);
   space_charge.set_fftw_helper(comm_spc, false); 
   space_charge.get_comm_sptr()->get_parent_sptr()->get();

   std::vector<Rectangular_grid_sptr>  Efield(space_charge.get_Efield(*rho_grid, bunch,3));
    
  
    double xnt=pi*nt/size[0];
    double ymt=pi*mt/size[1];
    double zpt=2.*pi*pt/size[2];
    double Cnmp=-(xnt*xnt+ymt*ymt+zpt*zpt);


    double sxnt=sin(pi*nt/double(shape[0]))/hi[0];
    double symt=sin(pi*mt/double(shape[1]))/hi[1];
    double szpt=sin(2.*pi*pt/double(shape[2]))/hi[2];


    MArray3d  rhocheck(boost::extents[shape[0]][shape[1]][shape[2]]);
  //  MArray3d  phi_check(boost::extents[shape[0]][shape[1]][shape[2]]);
    for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                for (int k = 0; k < shape[2]; ++k) {
                    rhocheck[i][j][k]=-rho[i][j][k]/(Cnmp*epsilon0);
                 //   phi_check[i][j][k] = phi_local->get_grid_points()[i][j][k]*phi_local->get_normalization();//epsilon0;
                }
            }
        }




 //  multi_array_check_equal(rhocheck, phi_check,  10.e-8);
  
    if (space_charge.get_comm_sptr()->has_this_rank()) {
	Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho_grid));
    
	Rectangular_grid_sptr  Ex(space_charge.get_En(*phi_local, 0));
	Rectangular_grid_sptr  Ey(space_charge.get_En(*phi_local, 1));
	Rectangular_grid_sptr  Ez(space_charge.get_En(*phi_local, 2));
	Ex->set_normalization(phi_local->get_normalization());
	Ey->set_normalization(phi_local->get_normalization());
	Ez->set_normalization(phi_local->get_normalization());
	BOOST_CHECK_EQUAL(Ex->get_domain().get_grid_shape()[0], space_charge.get_domain().get_grid_shape()[0]);
	BOOST_CHECK_EQUAL(Ex->get_domain().get_grid_shape()[1], space_charge.get_domain().get_grid_shape()[1]);
	BOOST_CHECK_EQUAL(Ex->get_domain().get_grid_shape()[2], space_charge.get_domain().get_grid_shape()[2]);
    
  
	MArray3d_ref Exa(Ex->get_grid_points());
	MArray3d_ref Eya(Ey->get_grid_points());
	MArray3d_ref Eza(Ez->get_grid_points());         
      
      
	for (int i = 0; i < shape[0]; ++i) {
		for (int j = 0; j < shape[1]; ++j) {
		    for (int k = 0; k < shape[2]; ++k) {
			Exa[i][j][k] *= Ex->get_normalization();
			Eya[i][j][k] *= Ey->get_normalization();
			Eza[i][j][k] *= Ez->get_normalization();			
		    }
		}
	    }
      // Check if Ex= -sxnt*cos(nt*(i+0.5)*pi/double(rho.shape()[0]))
      //       *sin(mt*(j+0.5)*pi/double(rho.shape()[1])) *cos((2.*pi*pt*k)/double(rho.shape()[2]))=
      //= -sxnt*cos(nt*(i+0.5)*pi/double(rho.shape()[0]))/ sin(nt*(i+0.5)*pi/double(rho.shape()[0]))*rho


      MArray3d  Excheck(boost::extents[shape[0]][shape[1]][shape[2]]);
      MArray3d  Eycheck(boost::extents[shape[0]][shape[1]][shape[2]]);
      MArray3d  Ezcheck(boost::extents[shape[0]][shape[1]][shape[2]]);
     
      for (int i = 0; i < shape[0]; ++i) {
		for (int j = 0; j < shape[1]; ++j) {
		    for (int k = 0; k < shape[2]; ++k) {
			Excheck[i][j][k]=Exa[i][j][k]*sin(nt*(i+0.5)*pi/double(shape[0]))/((-sxnt)*cos(nt*(i+0.5)*pi/double(shape[0])));
			Eycheck[i][j][k]=Eya[i][j][k]*sin(mt*(j+0.5)*pi/double(shape[1]))/((-symt)*cos(mt*(j+0.5)*pi/double(shape[1])));
			((k==shape[2]/2) || (k==0)) ? Ezcheck[i][j][k]= rhocheck[i][j][k]:
					  Ezcheck[i][j][k]=Eza[i][j][k]*cos((2.*pi*pt*k)/double(shape[2]))
					    /(szpt*sin((2.*pi*pt*k)/double(shape[2])));					    			

		      if ((i==0) || (i==shape[0]-1) ||(j==0) || (j==shape[1]-1)){
			Excheck[i][j][k]=0.;
			Eycheck[i][j][k]=0.;
			Ezcheck[i][j][k]=0.;		
			rhocheck[i][j][k]=0.;
		      }
		    }
		}
	    }
	if (space_charge.get_comm_sptr()->get_rank()==0) {   
	    multi_array_check_equal(rhocheck, Excheck,  10.e-8);
	    multi_array_check_equal(rhocheck, Eycheck,  10.e-8);
	    multi_array_check_equal(rhocheck, Ezcheck,  10.e-8);
	    
	}
    }//space_charge.get_comm_sptr()->has_this_rank()
     
    BOOST_CHECK_EQUAL(Efield[0]->get_domain().get_grid_shape()[0], space_charge.get_domain().get_grid_shape()[0]);
    BOOST_CHECK_EQUAL(Efield[0]->get_domain().get_grid_shape()[1], space_charge.get_domain().get_grid_shape()[1]);
    BOOST_CHECK_EQUAL(Efield[0]->get_domain().get_grid_shape()[2], space_charge.get_domain().get_grid_shape()[2]);
    BOOST_CHECK_EQUAL(Efield[1]->get_domain().get_grid_shape()[0], space_charge.get_domain().get_grid_shape()[0]);
    BOOST_CHECK_EQUAL(Efield[1]->get_domain().get_grid_shape()[1], space_charge.get_domain().get_grid_shape()[1]);
    BOOST_CHECK_EQUAL(Efield[1]->get_domain().get_grid_shape()[2], space_charge.get_domain().get_grid_shape()[2]);
    BOOST_CHECK_EQUAL(Efield[2]->get_domain().get_grid_shape()[0], space_charge.get_domain().get_grid_shape()[0]);
    BOOST_CHECK_EQUAL(Efield[2]->get_domain().get_grid_shape()[1], space_charge.get_domain().get_grid_shape()[1]);
    BOOST_CHECK_EQUAL(Efield[2]->get_domain().get_grid_shape()[2], space_charge.get_domain().get_grid_shape()[2]);
    MArray3d_ref Efieldx(Efield[0]->get_grid_points());
    MArray3d_ref Efieldy(Efield[1]->get_grid_points());
    MArray3d_ref Efieldz(Efield[2]->get_grid_points()); 
    
    for (int i = 0; i < shape[0]; ++i) {
	    for (int j = 0; j < shape[1]; ++j) {
		for (int k = 0; k < shape[2]; ++k) {		  
		    Efieldx[i][j][k] *= Efield[0]->get_normalization();
		    Efieldy[i][j][k] *= Efield[1]->get_normalization();
		    Efieldz[i][j][k] *= Efield[2]->get_normalization();
		}
	    }
    }
    MArray3d  Efxcheck(boost::extents[shape[0]][shape[1]][shape[2]]);
    MArray3d  Efycheck(boost::extents[shape[0]][shape[1]][shape[2]]);
    MArray3d  Efzcheck(boost::extents[shape[0]][shape[1]][shape[2]]);
    for (int i = 0; i < shape[0]; ++i) {
	for (int j = 0; j < shape[1]; ++j) {
	      for (int k = 0; k < shape[2]; ++k) {			
		  Efxcheck[i][j][k]=Efieldx[i][j][k]*sin(nt*(i+0.5)*pi/double(shape[0]))/((-sxnt)*cos(nt*(i+0.5)*pi/double(shape[0])));
		  Efycheck[i][j][k]=Efieldy[i][j][k]*sin(mt*(j+0.5)*pi/double(shape[1]))/((-symt)*cos(mt*(j+0.5)*pi/double(shape[1])));
		  ((k==shape[2]/2) || (k==0)) ? Efzcheck[i][j][k]= rhocheck[i][j][k]:
				    Efzcheck[i][j][k]=Efieldz[i][j][k]*cos((2.*pi*pt*k)/double(shape[2]))
				      /(szpt*sin((2.*pi*pt*k)/double(shape[2])));		
				      

		if ((i==0) || (i==shape[0]-1) ||(j==0) || (j==shape[1]-1)){			
		  Efxcheck[i][j][k]=0.;
		  Efycheck[i][j][k]=0.;
		  Efzcheck[i][j][k]=0.;
		  rhocheck[i][j][k]=0.;
		}
	      }
	  }
      } 
     
     multi_array_check_equal(rhocheck, Efxcheck,  10.e-8);
     multi_array_check_equal(rhocheck, Efycheck,  10.e-8);
     multi_array_check_equal(rhocheck, Efzcheck,  10.e-8);
}
     
     
BOOST_FIXTURE_TEST_CASE(get_En1, Ellipsoidal_bunch_fixture)
{

    std::vector<int > grid_shape1(3);
    grid_shape1[0] =32;
    grid_shape1[1] = 12;
    grid_shape1[2] = 10;
    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.1;
    size[2]=0.1;

    Space_charge_rectangular space_charge(size, grid_shape1);

    int nt=6;
    int mt=3;
    int pt=4;


    Rectangular_grid_sptr rho_grid(new Rectangular_grid(space_charge.get_domain_sptr()));
    std::vector<int > shape(space_charge.get_domain().get_grid_shape());

    MArray3d_ref rho(rho_grid->get_grid_points());
    for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                for (int k = 0; k < shape[2]; ++k) {
                   rho[i][j][k] = sin(nt*(i+0.5)*pi/double(shape[0]))*sin(mt*(j+0.5)*pi/double(shape[1])) *cos((2.*pi*pt*k)/double(shape[2]));
                }
            }
        }

   

   std::vector<double > hi(space_charge.get_domain().get_cell_size());
    int optimal_number=2;
    if (bunch.get_comm_sptr()->get_size() %  optimal_number==0) {
	  Commxx_divider comm_divider(optimal_number, false);
	  Commxx_sptr comm_spc=comm_divider.get_commxx_sptr(bunch.get_comm_sptr());
	  space_charge.set_fftw_helper(comm_spc, true);
	
	
	  double xnt=pi*nt/size[0];
	  double ymt=pi*mt/size[1];
	  double zpt=2.*pi*pt/size[2];
	  double Cnmp=-(xnt*xnt+ymt*ymt+zpt*zpt);


	  double sxnt=sin(pi*nt/double(shape[0]))/hi[0];
	  double symt=sin(pi*mt/double(shape[1]))/hi[1];
	  double szpt=sin(2.*pi*pt/double(shape[2]))/hi[2];


	  MArray3d  rhocheck(boost::extents[shape[0]][shape[1]][shape[2]]);
	//  MArray3d  phi_check(boost::extents[shape[0]][shape[1]][shape[2]]);
	  for (int i = 0; i < shape[0]; ++i) {
	      for (int j = 0; j < shape[1]; ++j) {
		  for (int k = 0; k < shape[2]; ++k) {
		      rhocheck[i][j][k]=-rho[i][j][k]/(Cnmp*epsilon0);
		  //   phi_check[i][j][k] = phi_local->get_grid_points()[i][j][k]*phi_local->get_normalization();//epsilon0;
		  }
	      }
	  }




	//  multi_array_check_equal(rhocheck, phi_check,  10.e-8);
    
	  Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho_grid));
      
	  Rectangular_grid_sptr  Ex(space_charge.get_En(*phi_local, 0));
	  Rectangular_grid_sptr  Ey(space_charge.get_En(*phi_local, 1));
	  Rectangular_grid_sptr  Ez(space_charge.get_En(*phi_local, 2));
	  BOOST_CHECK_EQUAL(Ex->get_domain().get_grid_shape()[0], space_charge.get_domain().get_grid_shape()[0]);
	  BOOST_CHECK_EQUAL(Ex->get_domain().get_grid_shape()[1], space_charge.get_domain().get_grid_shape()[1]);
	  BOOST_CHECK_EQUAL(Ex->get_domain().get_grid_shape()[2], space_charge.get_domain().get_grid_shape()[2]);
      
    
	  MArray3d_ref Exa(Ex->get_grid_points());
	  MArray3d_ref Eya(Ey->get_grid_points());
	  MArray3d_ref Eza(Ez->get_grid_points());         
	
	
	  for (int i = 0; i < shape[0]; ++i) {
		  for (int j = 0; j < shape[1]; ++j) {
		      for (int k = 0; k < shape[2]; ++k) {
			  Exa[i][j][k] *= Ex->get_normalization();
			  Eya[i][j][k] *= Ey->get_normalization();
			  Eza[i][j][k] *= Ez->get_normalization();			
		      }
		  }
	      }
	// Check if Ex= -sxnt*cos(nt*(i+0.5)*pi/double(rho.shape()[0]))
	//       *sin(mt*(j+0.5)*pi/double(rho.shape()[1])) *cos((2.*pi*pt*k)/double(rho.shape()[2]))=
	//= -sxnt*cos(nt*(i+0.5)*pi/double(rho.shape()[0]))/ sin(nt*(i+0.5)*pi/double(rho.shape()[0]))*rho


	  MArray3d  Excheck(boost::extents[shape[0]][shape[1]][shape[2]]);
	  MArray3d  Eycheck(boost::extents[shape[0]][shape[1]][shape[2]]);
	  MArray3d  Ezcheck(boost::extents[shape[0]][shape[1]][shape[2]]);
	
	  for (int i = 0; i < shape[0]; ++i) {
		    for (int j = 0; j < shape[1]; ++j) {
			for (int k = 0; k < shape[2]; ++k) {
			    Excheck[i][j][k]=Exa[i][j][k]*sin(nt*(i+0.5)*pi/double(shape[0]))/((-sxnt)*cos(nt*(i+0.5)*pi/double(shape[0])));
			    Eycheck[i][j][k]=Eya[i][j][k]*sin(mt*(j+0.5)*pi/double(shape[1]))/((-symt)*cos(mt*(j+0.5)*pi/double(shape[1])));
			    ((k==shape[2]/2) || (k==0)) ? Ezcheck[i][j][k]= rhocheck[i][j][k]:
					      Ezcheck[i][j][k]=Eza[i][j][k]*cos((2.*pi*pt*k)/double(shape[2]))
						/(szpt*sin((2.*pi*pt*k)/double(shape[2])));					    			

			  if ((i==0) || (i==shape[0]-1) ||(j==0) || (j==shape[1]-1)){
			    Excheck[i][j][k]=0.;
			    Eycheck[i][j][k]=0.;
			    Ezcheck[i][j][k]=0.;		
			    rhocheck[i][j][k]=0.;
			  }
			}
		    }
		}
	    
	  multi_array_check_equal(rhocheck, Excheck,  10.e-8);
	  multi_array_check_equal(rhocheck, Eycheck,  10.e-8);
	  multi_array_check_equal(rhocheck, Ezcheck,  10.e-8);
    }//	bunch.get_comm_sptr()->get_size() %  optimal_number    
     
}
   
BOOST_FIXTURE_TEST_CASE(serialize_, Ellipsoidal_bunch_fixture)
{
    
    std::vector<int > grid_shape1(3);
    grid_shape1[0] =32;
    grid_shape1[1] = 12;
    grid_shape1[2] = 10;
    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.1;
    size[2]=0.1;

   Space_charge_rectangular space_charge(size, grid_shape1);
   int optimal_number=3;
   Commxx_sptr comm_spc=make_optimal_spc_comm(bunch.get_comm_sptr(), optimal_number);
   space_charge.set_fftw_helper(comm_spc, false); 
    
 
   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string ss("space_charge");
    std::stringstream pp;
    pp<<rank;
    ss.append(pp.str());
    ss.append(".xml");
    xml_save(space_charge, ss);
  
    Space_charge_rectangular spc_loaded;
    xml_load(spc_loaded, ss);
 
    BOOST_CHECK_EQUAL(space_charge.get_pipe_size()[0], spc_loaded.get_pipe_size()[0]);
    BOOST_CHECK_EQUAL(space_charge.get_pipe_size()[1], spc_loaded.get_pipe_size()[1]);
    BOOST_CHECK_EQUAL(space_charge.get_pipe_size()[2], spc_loaded.get_pipe_size()[2]);
   
    BOOST_CHECK_EQUAL(space_charge.get_grid_shape()[0], spc_loaded.get_grid_shape()[0]);
    BOOST_CHECK_EQUAL(space_charge.get_grid_shape()[1], spc_loaded.get_grid_shape()[1]);
    BOOST_CHECK_EQUAL(space_charge.get_grid_shape()[2], spc_loaded.get_grid_shape()[2]);   
  // BOOST_CHECK_EQUAL(space_charge.have_fftw_helper,  spc_loaded.have_fftw_helper);
  // BOOST_CHECK_EQUAL(space_charge.use_comm_divider, spc_loaded.use_comm_divider);
    BOOST_CHECK_EQUAL(space_charge.get_comm_sptr()->has_this_rank(),spc_loaded.get_comm_sptr()->has_this_rank());
//    if (rank==0){
//       std::cout<<" MPI_IDENT="<<MPI_IDENT<<std::endl;
//       std::cout<<" MPI_CONGRUENT="<<MPI_CONGRUENT<<std::endl; 
//       std::cout<<" MPI_SIMILAR="<<MPI_SIMILAR<<std::endl;
//       std::cout<<" MPI_UNEQUAL="<<MPI_UNEQUAL<<std::endl;
//    }
    if (space_charge.get_comm_sptr()->has_this_rank()){
	int result;
	MPI_Comm_compare(space_charge.get_comm_sptr()->get(), spc_loaded.get_comm_sptr()->get(),  &result);
    //std::cout<<" result="<<result<<" rank="<<rank<<std::endl;  
	BOOST_CHECK((result == MPI_IDENT) || (result == MPI_CONGRUENT));

	int result_helper;
	MPI_Comm_compare(space_charge.get_fftw_helper_sptr()->get_comm_sptr()->get(), 
			     spc_loaded.get_fftw_helper_sptr()->get_comm_sptr()->get(),  &result_helper);
	BOOST_CHECK((result_helper == MPI_IDENT) || (result_helper == MPI_CONGRUENT));		     
	
    }
    
    
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_cell_size()[0],spc_loaded.get_domain_sptr()->get_cell_size()[0]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_cell_size()[1],spc_loaded.get_domain_sptr()->get_cell_size()[1]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_cell_size()[2],spc_loaded.get_domain_sptr()->get_cell_size()[2]);
   
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_physical_size()[0],spc_loaded.get_domain_sptr()->get_physical_size()[0]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_physical_size()[1],spc_loaded.get_domain_sptr()->get_physical_size()[1]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_physical_size()[2],spc_loaded.get_domain_sptr()->get_physical_size()[2]);
   
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_left()[0],spc_loaded.get_domain_sptr()->get_left()[0]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_left()[1],spc_loaded.get_domain_sptr()->get_left()[1]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_left()[2],spc_loaded.get_domain_sptr()->get_left()[2]);
   BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->is_periodic(), spc_loaded.get_domain_sptr()->is_periodic());  
    
}
 
 
BOOST_FIXTURE_TEST_CASE(serialize_comm_divider, Ellipsoidal_bunch_fixture)
{
    
    std::vector<int > grid_shape1(3);
    grid_shape1[0] =32;
    grid_shape1[1] = 12;
    grid_shape1[2] = 10;
    std::vector<double > size(3);
    size[0]=0.1;
    size[1]=0.1;
    size[2]=0.1;

   Space_charge_rectangular space_charge0(size, grid_shape1);
   int optimal_number=2;
   if (bunch.get_comm_sptr()->get_size() %  optimal_number==0) {
	Commxx_divider comm_divider(optimal_number, false);
	Commxx_sptr comm_spc=comm_divider.get_commxx_sptr(bunch.get_comm_sptr());
	space_charge0.set_fftw_helper(comm_spc, true);
        
	Space_charge_rectangular space_charge(*space_charge0.clone());
      
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::string ss("space_charge");
	std::stringstream pp;
	pp<<rank;
	ss.append(pp.str());
	ss.append(".xml");
	xml_save(space_charge, ss);
      
	Space_charge_rectangular spc_loaded;
	xml_load(spc_loaded, ss);
    
	BOOST_CHECK_EQUAL(space_charge.get_pipe_size()[0], spc_loaded.get_pipe_size()[0]);
	BOOST_CHECK_EQUAL(space_charge.get_pipe_size()[1], spc_loaded.get_pipe_size()[1]);
	BOOST_CHECK_EQUAL(space_charge.get_pipe_size()[2], spc_loaded.get_pipe_size()[2]);
      
	BOOST_CHECK_EQUAL(space_charge.get_grid_shape()[0], spc_loaded.get_grid_shape()[0]);
	BOOST_CHECK_EQUAL(space_charge.get_grid_shape()[1], spc_loaded.get_grid_shape()[1]);
	BOOST_CHECK_EQUAL(space_charge.get_grid_shape()[2], spc_loaded.get_grid_shape()[2]);   
	// BOOST_CHECK_EQUAL(space_charge.have_fftw_helper,  spc_loaded.have_fftw_helper);
	// BOOST_CHECK_EQUAL(space_charge.use_comm_divider, spc_loaded.use_comm_divider);
	BOOST_CHECK_EQUAL(space_charge.get_comm_sptr()->has_this_rank(),spc_loaded.get_comm_sptr()->has_this_rank());
	int result;
	MPI_Comm_compare(space_charge.get_comm_sptr()->get(), spc_loaded.get_comm_sptr()->get(),  &result);
	//std::cout<<" result="<<result<<" rank="<<rank<<std::endl;  
	BOOST_CHECK((result == MPI_IDENT) || (result == MPI_CONGRUENT));  
	
	int result_helper;
	MPI_Comm_compare(space_charge.get_fftw_helper_sptr()->get_comm_sptr()->get(), 
			     spc_loaded.get_fftw_helper_sptr()->get_comm_sptr()->get(),  &result_helper);
	BOOST_CHECK((result_helper == MPI_IDENT) || (result_helper == MPI_CONGRUENT));		
	
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_cell_size()[0],spc_loaded.get_domain_sptr()->get_cell_size()[0]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_cell_size()[1],spc_loaded.get_domain_sptr()->get_cell_size()[1]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_cell_size()[2],spc_loaded.get_domain_sptr()->get_cell_size()[2]);
	
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_physical_size()[0],spc_loaded.get_domain_sptr()->get_physical_size()[0]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_physical_size()[1],spc_loaded.get_domain_sptr()->get_physical_size()[1]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_physical_size()[2],spc_loaded.get_domain_sptr()->get_physical_size()[2]);
	
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_left()[0],spc_loaded.get_domain_sptr()->get_left()[0]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_left()[1],spc_loaded.get_domain_sptr()->get_left()[1]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_left()[2],spc_loaded.get_domain_sptr()->get_left()[2]);
	BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->is_periodic(), spc_loaded.get_domain_sptr()->is_periodic()); 
	
   }			
}
 
