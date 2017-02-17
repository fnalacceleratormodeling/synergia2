#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/diagnostics_apertures_loss.h"
#include "synergia/lattice/lattice.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/period.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
const char attribute[] = "test";
double const tolerance1=1.e-16;
struct Fixture
{
    Fixture() :
        lattice_sptr(new Lattice("lattice")), f("quadrupole", "f"),
                o("drift", "o"), d("quadrupole", "d")
    {
        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        BOOST_TEST_MESSAGE("setup fixture");
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Lattice_sptr lattice_sptr;
    Lattice_element f, o, d;
};

const double real_num = 1.7e11;
const int total_num = 10000;
const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double total_energy = 1.1;

struct Fixture_bunch
{
    Fixture_bunch() :
        lattice_sptr(new Lattice("lattice")), f("quadrupole", "f"),
                o("drift", "o"), d("quadrupole", "d"), comm_sptr(new Commxx), four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum),bunch(reference_particle,  total_num, real_num, comm_sptr),
                seed(718281828), distribution(seed, *comm_sptr)
              
    {
        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        stdx = 0.2;
        stdy = 0.03;
        stdz = 1.;
        covariances[0][0] = stdx * stdx;
        covariances[2][2] = stdy * stdy;
        covariances[4][4] = stdz * stdz;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 0.00001;
        populate_6d(distribution, bunch, means, covariances);
        bunch.set_bucket_index(12);
        BOOST_TEST_MESSAGE("setup bunch_lattice_fixture");
    }


    
    
    
    ~Fixture_bunch()
    {
        BOOST_TEST_MESSAGE("teardown bunch_lattice_fixture");
    }

    Lattice_sptr lattice_sptr;
    Lattice_element f, o, d;
    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    Random_distribution distribution;
    double stdx, stdy, stdz;
};

BOOST_AUTO_TEST_CASE(construct_diagnostics)
{

   Diagnostics_loss_sptr diag_loss_sptr(new Diagnostics_loss("aperture_loss_file.h5","aperture"));
   BOOST_CHECK_EQUAL(diag_loss_sptr->get_type(),Diagnostics_loss::aperture_type);
   Diagnostics_loss_sptr diag_zcut_sptr(new Diagnostics_loss("zcut_loss_file.h5","zcut"));  
   BOOST_CHECK_EQUAL(diag_zcut_sptr->get_type(),Diagnostics_loss::zcut_type);

}

 
BOOST_FIXTURE_TEST_CASE(set_diagnostics, Fixture)
{
     Diagnostics_loss_sptr 
     diag_loss_sptr=Diagnostics_loss_sptr(new Diagnostics_loss("aperture_loss_file.h5","aperture"));
     lattice_sptr->add_loss_diagnostics(diag_loss_sptr);
   
}

BOOST_FIXTURE_TEST_CASE(update_diagnostics, Fixture_bunch)
{
     Diagnostics_loss_sptr 
     diag_loss_sptr=Diagnostics_loss_sptr(new Diagnostics_loss("aperture_loss_file.h5","aperture"));
     
     
     lattice_sptr->add_loss_diagnostics(diag_loss_sptr);
     
     Diagnostics_loss_sptr  diag_sptr;
     diag_sptr= *(lattice_sptr->get_loss_diagnostics_list().begin());
     
     Bunch_sptr bunch_sptr((new Bunch(bunch)));
     diag_sptr->set_bunch_sptr(bunch_sptr);
     Commxx_sptr comm_sptr(bunch_sptr->get_comm_sptr());
     

     int total_size=comm_sptr->get_size();
     BOOST_CHECK_EQUAL(total_size, 3);
     int rank=comm_sptr->get_rank();

     MArray1d coords(boost::extents[6]);
     int b_index=10;
     int repetition=3;
     double s=100.10;
     double s_n=5.4;
      coords[0]=1.2+10*rank;
      coords[1]=2.2+10*rank;
      coords[2]=3.2+10*rank;
      coords[3]=4.2+10*rank;
      coords[4]=5.2+10*rank;
      coords[5]=6.2+10*rank;                            
      diag_sptr->update( b_index, repetition, s,s_n, coords );
      if (rank==0){       
          b_index=20;
          repetition=6;
          s=200.10;
          s_n=8.4;
          coords[0]=-1.2;
          coords[1]=-2.2;
          coords[2]=-3.2;
          coords[3]=-4.2;
          coords[4]=-5.2;
          coords[5]=-6.2; 
          diag_sptr->update( b_index, repetition, s,s_n, coords );   
      }
      
      diag_sptr->write();
      
      if (rank==1){
          b_index=30;
          repetition=8;
          s=300.10;
          s_n=12.4;
          coords[0]=-10.2;
          coords[1]=-20.2;
          coords[2]=-30.2;
          coords[3]=-40.2;
          coords[4]=-50.2;
          coords[5]=-60.2; 
          diag_sptr->update( b_index, repetition, s,s_n, coords );   
      }
      
      diag_sptr->write();
      
      if ( diag_sptr->get_write_helper().write_locally()) {
         diag_sptr->get_write_helper().get_hdf5_file_sptr()->close();
         Hdf5_file file("aperture_loss_file.h5", Hdf5_file::read_only); 
         MArray1i read_rep=file.read<MArray1i > ("repetition");
         MArray1i read_bi=file.read<MArray1i > ("bucket_index");
         MArray1d read_s=file.read<MArray1d > ("s");
         MArray1d read_sn=file.read<MArray1d > ("s_n");
         MArray2d read_cc=file.read<MArray2d > ("coordinates");
         file.close();
         BOOST_CHECK_EQUAL(read_rep.size(), 5);
         BOOST_CHECK_EQUAL(read_bi.size(), 5);
         BOOST_CHECK_EQUAL(read_s.size(), 5);
         BOOST_CHECK_EQUAL(read_sn.size(), 5);
         BOOST_CHECK_EQUAL(read_cc.shape()[0], 6);
         BOOST_CHECK_EQUAL(read_cc.shape()[1], 5);
         BOOST_CHECK_CLOSE(read_cc[0][0],1.2, tolerance1 );
         BOOST_CHECK_CLOSE(read_cc[1][0],2.2, tolerance1 );
         BOOST_CHECK_CLOSE(read_cc[2][0],3.2, tolerance1 );
         BOOST_CHECK_CLOSE(read_cc[5][0],6.2, tolerance1 ); 
         BOOST_CHECK_CLOSE(read_cc[4][1],-5.2, tolerance1 ); 
         BOOST_CHECK_CLOSE(read_cc[0][4],-10.2, tolerance1 ); 
         BOOST_CHECK_CLOSE(read_s[0],100.10, tolerance1 ); 
         BOOST_CHECK_CLOSE(read_sn[1],8.4, tolerance1 ); 
         BOOST_CHECK_EQUAL(read_bi[4],30); 
         BOOST_CHECK_EQUAL(read_rep[4],8); 
         
      }
      
      
}

BOOST_FIXTURE_TEST_CASE(zcut, Fixture_bunch)
{    
      
    Diagnostics_loss_sptr 
     diag_loss_sptr=Diagnostics_loss_sptr(new Diagnostics_loss("zcut_loss_file.h5","zcut"));
    // lattice_sptr->add_apertures_diagnostics(diag_loss_sptr);
     Bunch_sptr bunch_sptr((new Bunch(bunch)));
     int ini_bch_num=bunch_sptr->get_total_num();
     diag_loss_sptr->set_bunch_sptr(bunch_sptr);
          
     
     double length=2*stdz;
     double beta=bunch_sptr->get_reference_particle().get_beta();
     apply_zcut(*bunch_sptr, length,diag_loss_sptr );
     diag_loss_sptr->write();
   
     
      MArray2d_ref particles(bunch_sptr->get_local_particles());
     
     for (int p=0;p<bunch_sptr->get_local_num();++p){        
         BOOST_CHECK(fabs(particles[p][4])<0.5*length/beta);
      }
     
     
     
      if ( diag_loss_sptr->get_write_helper().write_locally()) {
         diag_loss_sptr->get_write_helper().get_hdf5_file_sptr()->close();
          Hdf5_file file("zcut_loss_file.h5", Hdf5_file::read_only); 
          MArray1i read_rep=file.read<MArray1i > ("repetition");
          MArray1i read_bi=file.read<MArray1i > ("bucket_index");
          MArray1d read_s=file.read<MArray1d > ("s");
          MArray1d read_sn=file.read<MArray1d > ("s_n");
          MArray2d read_cc=file.read<MArray2d > ("coordinates");
          file.close();
         double nloss=read_rep.size();
         BOOST_CHECK_EQUAL(nloss,ini_bch_num-bunch_sptr->get_total_num());
         for (int p=0;p<nloss;++p){
            BOOST_CHECK_EQUAL(read_bi[p],12);
            BOOST_CHECK(fabs(read_cc[4][p])>=0.5*length/beta);
         }             
      }     
}