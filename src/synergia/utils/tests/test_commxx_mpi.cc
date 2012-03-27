#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/commxx.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct1)
{
    Commxx commxx;
}

BOOST_AUTO_TEST_CASE(construct2)
{
    bool per_host = true;
    Commxx commxx(per_host);
}

BOOST_AUTO_TEST_CASE(construct3)
{
    Commxx_sptr parent_sptr(new Commxx);
    std::vector<int > ranks(1);
    ranks.at(0) = 0;
    Commxx commxx(parent_sptr, ranks);
}

BOOST_AUTO_TEST_CASE(construct4)
{
    Commxx_sptr parent_sptr(new Commxx);
    std::vector<int > ranks(1);
    ranks.at(0) = 0;
    bool per_host = true;
    Commxx(parent_sptr, ranks, per_host);
}

BOOST_AUTO_TEST_CASE(get_rank)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int commxx_rank = Commxx().get_rank();
    BOOST_CHECK_EQUAL(rank, commxx_rank);
}

BOOST_AUTO_TEST_CASE(get_size)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int commxx_size = Commxx().get_size();
    BOOST_CHECK_EQUAL(size, commxx_size);
}

BOOST_AUTO_TEST_CASE(has_this_rank1)
{
    Commxx commxx;
    BOOST_CHECK(commxx.has_this_rank());
}

BOOST_AUTO_TEST_CASE(has_this_rank2)
{
    Commxx commxx(true);
    BOOST_CHECK(commxx.has_this_rank());
}

BOOST_AUTO_TEST_CASE(has_this_rank3)
{
    Commxx_sptr parent_sptr(new Commxx);
    std::vector<int > ranks(1);
    const int only_rank = 0;
    ranks.at(0) = only_rank;
    Commxx commxx(parent_sptr, ranks);

    BOOST_CHECK_EQUAL(commxx.has_this_rank(), (parent_sptr->get_rank() == only_rank));
}

BOOST_AUTO_TEST_CASE(get)
{
    MPI_Comm comm = Commxx().get();
    int result;
    MPI_Comm_compare(MPI_COMM_WORLD, comm, &result);
    BOOST_CHECK(result == MPI_IDENT);
}

BOOST_AUTO_TEST_CASE(serialize1)
{
    const std::string serialize_file_name("commxx1.xml");
    {
        Commxx commxx;
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

BOOST_AUTO_TEST_CASE(serialize2)
{
    const std::string serialize_file_name("commxx2.xml");
    {
        Commxx commxx(true);
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

BOOST_AUTO_TEST_CASE(serialize3)
{
    const std::string serialize_file_name("commxx3.xml");
    {
        Commxx_sptr parent_sptr(new Commxx);
        std::vector<int > ranks(1);
        ranks.at(0) = 0;
        Commxx commxx(parent_sptr, ranks);
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

BOOST_AUTO_TEST_CASE(serialize4)
{
    const std::string serialize_file_name("commxx4.xml");
    {
        Commxx_sptr parent_sptr(new Commxx);
        std::vector<int > ranks(1);
        ranks.at(0) = 0;
        bool per_host = true;
        Commxx commxx(parent_sptr, ranks, per_host);
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}
