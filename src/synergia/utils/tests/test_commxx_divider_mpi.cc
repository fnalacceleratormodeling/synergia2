#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/commxx_divider.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

BOOST_AUTO_TEST_CASE(construct1)
{
    Commxx_divider commxx_divider;
}

BOOST_AUTO_TEST_CASE(construct2)
{
    int subsize = 3;
    bool per_host = true;
    Commxx_divider commxx_divider(subsize, per_host);
}

BOOST_AUTO_TEST_CASE(get_commxx_sptr1)
{
    int subsize = 1;
    bool per_host = false;
    Commxx_divider commxx_divider(subsize, per_host);
    Commxx_sptr parent(new Commxx());
    BOOST_CHECK(commxx_divider.get_commxx_sptr(parent)->get_size() == subsize);
}

BOOST_AUTO_TEST_CASE(get_commxx_sptr2)
{
    int subsize = 2;
    bool per_host = false;
    Commxx_divider commxx_divider(subsize, per_host);
    Commxx_sptr parent(new Commxx());
    if (subsize < parent->get_size()) {
        if (parent->get_size() % subsize == 0) {
            BOOST_CHECK(
                    commxx_divider.get_commxx_sptr(parent)->get_size() == subsize);
        } else {
            bool caught = false;
            try {
                Commxx_sptr result = commxx_divider.get_commxx_sptr(parent);
            }
            catch (std::runtime_error &) {
                caught = true;
            }
            BOOST_CHECK(caught);
        }
    } else {
        BOOST_CHECK(
                commxx_divider.get_commxx_sptr(parent)->get_size() == parent->get_size());
    }
}

BOOST_AUTO_TEST_CASE(get_commxx_sptr3)
{
    int subsize = 3;
    bool per_host = false;
    Commxx_divider commxx_divider(subsize, per_host);
    Commxx_sptr parent(new Commxx());
    if (subsize < parent->get_size()) {
        if (parent->get_size() % subsize == 0) {
            BOOST_CHECK(
                    commxx_divider.get_commxx_sptr(parent)->get_size() == subsize);
        } else {
            bool caught = false;
            try {
                Commxx_sptr result = commxx_divider.get_commxx_sptr(parent);
            }
            catch (std::runtime_error &) {
                caught = true;
            }
            BOOST_CHECK(caught);
        }
    } else {
        BOOST_CHECK(
                commxx_divider.get_commxx_sptr(parent)->get_size() == parent->get_size());
    }
}

BOOST_AUTO_TEST_CASE(get_commxx_sptr4)
{
    int subsize = 4;
    bool per_host = false;
    Commxx_divider commxx_divider(subsize, per_host);
    Commxx_sptr parent(new Commxx());
    if (subsize < parent->get_size()) {
        if (parent->get_size() % subsize == 0) {
            BOOST_CHECK(
                    commxx_divider.get_commxx_sptr(parent)->get_size() == subsize);
        } else {
            bool caught = false;
            try {
                Commxx_sptr result = commxx_divider.get_commxx_sptr(parent);
            }
            catch (std::runtime_error &) {
                caught = true;
            }
            BOOST_CHECK(caught);
        }
    } else {
        BOOST_CHECK(
                commxx_divider.get_commxx_sptr(parent)->get_size() == parent->get_size());
    }
}

BOOST_AUTO_TEST_CASE(get_commxx_sptr_cached)
{
    int subsize = 1;
    bool per_host = false;
    Commxx_divider commxx_divider(subsize, per_host);
    Commxx_sptr parent(new Commxx());
    Commxx_sptr initial(commxx_divider.get_commxx_sptr(parent));
    Commxx_sptr cached(commxx_divider.get_commxx_sptr(parent));
    BOOST_CHECK_EQUAL(initial, cached);
}

BOOST_AUTO_TEST_CASE(serialize_)
{
    const std::string serialize_file_name("commxx_divider.xml");
    {
        Commxx_divider commxx_divider(7, true);
        remove_serialization_directory();
        xml_save(commxx_divider,
                get_serialization_path(serialize_file_name).c_str(), true);
    }

    {
        Commxx_divider commxx_divider_resumed;
        xml_load(commxx_divider_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}
