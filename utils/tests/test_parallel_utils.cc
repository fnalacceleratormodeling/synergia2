#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/boost_test_mpi_fixture.h"
#include "utils/parallel_utils.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(test_Commxx_construct)
{
    Commxx comm(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_Commxx_get_rank)
{
    int rank = Commxx(MPI_COMM_WORLD).get_rank();
    BOOST_CHECK_EQUAL(rank,0);
}

BOOST_AUTO_TEST_CASE(test_Commxx_get_size)
{
    int size = Commxx(MPI_COMM_WORLD).get_size();
    BOOST_CHECK_EQUAL(size,1);
}

BOOST_AUTO_TEST_CASE(test_Commxx_get)
{
    MPI_Comm comm = Commxx(MPI_COMM_WORLD).get();
    BOOST_CHECK_EQUAL(comm,MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_decompose_1d_raw1)
{
    // Test counts and offsets when length is a multiple of
    // number of processors
    const int procs = 4;
    std::vector<int > offsets(procs), counts(procs);
    const int num_per_proc = 5;
    decompose_1d_raw(procs, procs * num_per_proc, offsets, counts);
    for (int i = 0; i < procs; ++i) {
        BOOST_CHECK_EQUAL(num_per_proc,counts.at(i));
        BOOST_CHECK_EQUAL(num_per_proc*i,offsets.at(i));
    }
}

BOOST_AUTO_TEST_CASE(test_decompose_1d_raw2)
{
    // Test counts and offsets when length is not a multiple of
    // number of processors
    const int procs = 4;
    std::vector<int > offsets(procs), counts(procs);
    const int num_per_proc = 5;
    const int remainder = 3;
    int total = procs * num_per_proc + remainder;
    decompose_1d_raw(procs, total, offsets, counts);
    int min = num_per_proc * 2;
    int max = 0;
    int expected_offset = 0;
    for (int i = 0; i < procs; ++i) {
        BOOST_CHECK_EQUAL(expected_offset,offsets.at(i));
        int count = counts.at(i);
        if (count < min) {
            min = count;
        }
        if (count > max) {
            max = count;
        }
        expected_offset += count;
    }
    // The zeroth processor should always get the minimum number
    BOOST_CHECK_EQUAL(min,counts.at(0));
    // No two processors should differ by more than one
    BOOST_CHECK((max-min) <= 1);
    // The expected_offset at the end is the sum
    BOOST_CHECK_EQUAL(expected_offset,total);
}

BOOST_AUTO_TEST_CASE(test_decompose_1d_raw3)
{
    // Test counts and offsets when length is less than
    // number of processors
    const int procs = 4;
    std::vector<int > offsets(procs), counts(procs);
    const int num_per_proc = 5;
    const int empty = 2;
    int total = procs * num_per_proc - empty;
    decompose_1d_raw(procs, total, offsets, counts);
    int min = num_per_proc * 2;
    int max = 0;
    int expected_offset = 0;
    for (int i = 0; i < procs; ++i) {
        BOOST_CHECK_EQUAL(expected_offset,offsets.at(i));
        int count = counts.at(i);
        if (count < min) {
            min = count;
        }
        if (count > max) {
            max = count;
        }
        expected_offset += count;
    }
    // The zeroth processor should always get the minimum number
    BOOST_CHECK_EQUAL(min,counts.at(0));
    // No two processors should differ by more than one
    BOOST_CHECK((max-min) <= 1);
    // The expected_offset at the end is the sum
    BOOST_CHECK_EQUAL(expected_offset,total);
}

BOOST_AUTO_TEST_CASE(test_decompose_1d)
{
    const int procs = 1;
    std::vector<int > offsets(procs), counts(procs);
    const int length = 17;
    decompose_1d(Commxx(MPI_COMM_WORLD), length, offsets, counts);
    BOOST_CHECK_EQUAL(0,offsets.at(0));
    BOOST_CHECK_EQUAL(length,counts.at(0));
}

BOOST_AUTO_TEST_CASE(test_decompose_1d_local)
{
    const int length = 17;
    int local_length =
            decompose_1d_local(Commxx(MPI_COMM_WORLD), length);
    BOOST_CHECK_EQUAL(length,local_length);
}
