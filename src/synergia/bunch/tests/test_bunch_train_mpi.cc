#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/train.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct1)
{
    const int num_bunches = 1;
    Commxx mcommx=Commxx(MPI_COMM_WORLD);
    Train_comms train(num_bunches, mcommx);
}

BOOST_AUTO_TEST_CASE(construct2)
{

    const int num_bunches = 2;
    Commxx mcommx=Commxx(MPI_COMM_WORLD);
    Train_comms train(num_bunches, mcommx);
}   

BOOST_AUTO_TEST_CASE(construct3)
{
    const int num_bunches = 16;
    Commxx mcommx=Commxx(MPI_COMM_WORLD);
    Train_comms train(num_bunches, mcommx);
   
}

BOOST_AUTO_TEST_CASE(get_master_comm)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
    BOOST_CHECK_EQUAL(bunch_train.get_master_comm().get(), MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(get_comm)
{
    const int num_bunches = 16;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        bunch_train.get_comm(bunch);
        // don't know what to check here...
    }
}

BOOST_AUTO_TEST_CASE(is_on_this_rank)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        bool on_this_rank = bunch_train.is_on_this_rank(bunch);
        // could use a better test here...
    }
}

BOOST_AUTO_TEST_CASE(diag_get_comm)
{
    const int num_bunches = 2;
    const double bunch_separation = 1.7;
    Commxx mcommx=Commxx(MPI_COMM_WORLD);
//     std::cout<<" mcommx rank="<<mcommx.get_rank()<<std::endl;
//     std::cout<<" mcommx size="<<mcommx.get_size()<<std::endl;
    
    Train_comms train(num_bunches, mcommx);
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        if (train.is_on_this_rank(bunch)) {
           int  rk, sz;
            rk = train.get_comm(bunch).get_rank();
            sz= train.get_comm(bunch).get_size();
           // std::cout<<" bunch="<<bunch<<" comm rank="<<rk<<" size="<<sz<<" master com rank="<<mcommx.get_rank()<<std::endl;
       //  std::cout<<" comm size="<<comm.get_size()<<std::endl;
        }
     
     
    }
}
