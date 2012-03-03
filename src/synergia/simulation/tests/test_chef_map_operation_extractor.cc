#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/operation_extractor.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_map_operation_extractor chef_map_o_e(chef_lattice_sptr, map_order);
}

BOOST_FIXTURE_TEST_CASE(extract1, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_map_operation_extractor chef_map_o_e(chef_lattice_sptr, map_order);
    Lattice_elements elements(lattice_sptr->get_elements());
    double fraction = 0.5;
    double slice_end = elements.front()->get_length() * fraction;
    Lattice_element_slice_sptr slice_sptr(new Lattice_element_slice(
            elements.front(), 0, slice_end));
    Lattice_element_slices slices;
    slices.push_back(slice_sptr);
    chef_lattice_sptr->construct_sliced_beamline(slices);
    Independent_operations ops = chef_map_o_e.extract(b.reference_particle,
            slices);
    BOOST_CHECK_EQUAL(ops.size(), 1);
    BOOST_CHECK_EQUAL(ops.front()->get_type(), fast_mapping_type_name);
}
// test_note: we still need to check the extracted value

BOOST_FIXTURE_TEST_CASE(extract3, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_map_operation_extractor chef_map_o_e(chef_lattice_sptr, map_order);
    Lattice_elements elements(lattice_sptr->get_elements());

    Lattice_element_slices slices;
    Lattice_elements::iterator it = elements.begin();
    // get three slices...

    // ... the last first_fraction of the first element
    double first_fraction = 0.2;
    double first_length = (*it)->get_length();
    double first_begin = (1 - first_fraction) * first_length;
    Lattice_element_slice_sptr slice_sptr1(new Lattice_element_slice(*it,
            first_begin, first_length));
    slices.push_back(slice_sptr1);

    // ... the entire second element
    ++it;
    Lattice_element_slice_sptr slice_sptr2(new Lattice_element_slice(*it));
    slices.push_back(slice_sptr2);

    /// the first third_fraction of the third element
    ++it;
    double third_fraction = 0.3;
    double third_length = (*it)->get_length();
    double third_end = third_fraction * third_length;
    Lattice_element_slice_sptr slice_sptr3(new Lattice_element_slice(*it, 0,
            third_end));
    slices.push_back(slice_sptr3);

    chef_lattice_sptr->construct_sliced_beamline(slices);
    Independent_operations ops = chef_map_o_e.extract(b.reference_particle,
            slices);
    BOOST_CHECK_EQUAL(ops.size(), 1);
    BOOST_CHECK_EQUAL(ops.front()->get_type(), fast_mapping_type_name);
}
// test_note: we still need to check the extracted value

BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_map_operation_extractor chef_map_o_e(chef_lattice_sptr, map_order);
    xml_save(chef_map_o_e, "chef_map_operation_extractor.xml");

    Chef_map_operation_extractor loaded;
    xml_load(loaded, "chef_map_operation_extractor.xml");
}
