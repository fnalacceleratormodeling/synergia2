#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/operation_extractor.h"
#include "components/lattice/lattice_element_slice.h"
#include "lattice_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(*lattice_sptr));
    Chef_map_operation_extractor chef_map_o_e(chef_lattice_sptr, map_order);
}
//extract(Reference_particle const& reference_particle,
//        Lattice_element_slices const& slices);

BOOST_FIXTURE_TEST_CASE(extract1, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(*lattice_sptr));
    Chef_map_operation_extractor chef_map_o_e(chef_lattice_sptr, map_order);
    Lattice_elements elements(lattice_sptr->get_elements());
    double fraction = 0.5;
    double slice_end = elements.front().get_length()*fraction;
    Lattice_element_slice_sptr slice_sptr(new Lattice_element_slice(elements.front(),0,slice_end));
    Lattice_element_slices slices;
    slices.push_back(slice_sptr);
    std::cout << "jfa in test 1\n";
    chef_lattice_sptr->construct_sliced_beamline(slices);
    std::cout << "jfa in test 2\n";
    chef_map_o_e.extract(b.reference_particle,slices);
}

//BOOST_FIXTURE_TEST_CASE(get_type, Chef_elements_fixture)
//{
//    Chef_propagate_operation chef_propagate_operation(chef_elements);
//    BOOST_CHECK_EQUAL(chef_propagate_operation.get_type(), "chef_propagate");
//}
//
//
//BOOST_FIXTURE_TEST_CASE(apply, Chef_elements_fixture)
//{
//    Bunch_fixture b;
//
//    Chef_propagate_operation chef_propagate_operation(chef_elements);
//    //    multi_array_print(b.bunch.get_local_particles(), "particles before");
//    chef_propagate_operation.apply(b.bunch);
//    //    multi_array_print(b.bunch.get_local_particles(), "particles after");
//}
// test_note: We need to check that apply actual produces the correct results.
//            As of this writing, it almost certainly doesn't
