#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/parallel_utils.h"
#include <cmath>
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(test_Commxx_construct1)
{
    Commxx comm(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_Commxx_construct2)
{
    Commxx
    comm();
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

BOOST_AUTO_TEST_CASE(test_Commxx_set_get)
{
    MPI_Comm new_comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &new_comm);
    Commxx commxx;
    commxx.set(new_comm);
    MPI_Comm comm = commxx.get();
    BOOST_CHECK(comm != MPI_COMM_WORLD);
    BOOST_CHECK_EQUAL(comm, new_comm);
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
    int local_length = decompose_1d_local(Commxx(MPI_COMM_WORLD), length);
    BOOST_CHECK_EQUAL(length,local_length);
}

void
verify_ranks_elms_gt_procs(int elements, int procs,
        std::vector<std::vector<int > > const& ranks)
{
    BOOST_CHECK_EQUAL(ranks.size(), elements);
    std::vector<int > elements_on_proc(procs);
    for (int i = 0; i < procs; ++i) {
        elements_on_proc.at(i) = 0;
    }
    for (int element = 0; element < elements; ++element) {
        BOOST_CHECK_EQUAL((ranks.at(element)).size(), 1);
        ++elements_on_proc.at(ranks.at(element).at(0));
    }
    int max_on_proc = elements_on_proc.at(0);
    int min_on_proc = elements_on_proc.at(0);
    for (int i = 1; i < procs; ++i) {
        if (elements_on_proc.at(i) > max_on_proc) {
            max_on_proc = elements_on_proc.at(i);
        }
        if (elements_on_proc.at(i) < min_on_proc) {
            min_on_proc = elements_on_proc.at(i);
        }
    }
    BOOST_CHECK_EQUAL(max_on_proc, std::ceil((1.0*elements)/procs));
    BOOST_CHECK_EQUAL(min_on_proc, std::floor(elements/procs));
}

BOOST_AUTO_TEST_CASE(test_distribute_1d_raw1)
{
    const int elements = 4;
    const int procs = 1;
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    verify_ranks_elms_gt_procs(elements, procs, ranks);
}

BOOST_AUTO_TEST_CASE(test_distribute_1d_raw2)
{
    const int elements = 8;
    const int procs = 2;
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    verify_ranks_elms_gt_procs(elements, procs, ranks);
}

BOOST_AUTO_TEST_CASE(test_distribute_1d_raw3)
{
    const int elements = 57;
    const int procs = 13;
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    verify_ranks_elms_gt_procs(elements, procs, ranks);
}

void
verify_ranks_elms_lt_procs(int elements, int procs,
        std::vector<std::vector<int > > const& ranks)
{
    BOOST_CHECK_EQUAL(ranks.size(), elements);
    std::vector<bool > used(procs);
    for (int proc = 0; proc < procs; ++proc) {
        used.at(proc) = false;
    }
    int max_size = ranks.at(0).size();
    int min_size = ranks.at(0).size();
    for (int element = 0; element < elements; ++element) {
        if (ranks.at(element).size() > max_size) {
            max_size = ranks.at(element).size();
        }
        if (ranks.at(element).size() < max_size) {
            min_size = ranks.at(element).size();
        }
        for (int i = 0; i < ranks.at(element).size(); ++i) {
            BOOST_CHECK(!used.at(ranks.at(element).at(i)));
            used.at(ranks.at(element).at(i)) = true;
        }
    }
    BOOST_CHECK_EQUAL(max_size, std::ceil((1.0*procs)/elements));
    BOOST_CHECK_EQUAL(min_size, std::floor(procs/elements));
}

BOOST_AUTO_TEST_CASE(test_distribute_1d_raw4)
{
    const int elements = 1;
    const int procs = 32;
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    verify_ranks_elms_lt_procs(elements, procs, ranks);
}

BOOST_AUTO_TEST_CASE(test_distribute_1d_raw5)
{
    const int elements = 8;
    const int procs = 32;
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    verify_ranks_elms_lt_procs(elements, procs, ranks);
}

BOOST_AUTO_TEST_CASE(test_distribute_1d_raw6)
{
    const int elements = 17;
    const int procs = 71;
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    verify_ranks_elms_lt_procs(elements, procs, ranks);
}

BOOST_AUTO_TEST_CASE(test_distribute_1d)
{
    const int elements = 17;
    std::vector<std::vector<int > > ranks(
            distribute_1d(Commxx(MPI_COMM_WORLD), elements));
    verify_ranks_elms_gt_procs(elements, 1, ranks);
}
/*
BOOST_AUTO_TEST_CASE(test_print_distribute_1d_raw1)
{   
//     const int elements = 6;
//     const int procs = 32;
            
    const int elements = 5;
    const int procs = 12;       
            
    std::vector<std::vector<int > > ranks(distribute_1d_raw(procs, elements));
    BOOST_CHECK_EQUAL(ranks.size(), elements);        
    std::cout<<" ranks size="<<ranks.size()<<std::endl;
    for (int element = 0; element < elements; ++element) { 
        std::cout<<" element="<<element<<" ranks size at element="<<ranks.at(element).size()<<std::endl;
        std::cout<<"       ******** procs for element=";
        for (int j = 0; j < ranks.at(element).size(); ++j) { 
            std::cout<<ranks.at(element).at(j)<<"  ";
        }
        std::cout<<std::endl;
    }       
}*/

BOOST_AUTO_TEST_CASE(test_counts_and_offsets)
{   
//     const int elements = 5;
//     const int procs = 12;
            
    const int elements = 5;
    const int procs = 12;       
    std::vector<int > counts(procs);
    std::vector<int > offsets(procs);        
    counts_and_offsets_for_impedance_raw(procs, elements,  offsets, counts);
  //  for (int i = 0; i < procs; ++i){
 //       std::cout<<" proc="<<i<<" count="<<counts[i]<<" offset="<<offsets[i]<<std::endl;
 //   }
            
   
}

