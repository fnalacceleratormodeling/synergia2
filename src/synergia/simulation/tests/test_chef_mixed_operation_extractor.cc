#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/operation_extractor.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_mixed_operation_extractor mixed_chef_o_e(chef_lattice_sptr, map_order);
}

BOOST_FIXTURE_TEST_CASE(extract1, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_mixed_operation_extractor mixed_chef_o_e(chef_lattice_sptr, map_order);
    Lattice_elements elements(lattice_sptr->get_elements());
    double fraction = 0.5;
    double slice_end = elements.front()->get_length() * fraction;
    Lattice_element_slice_sptr slice_sptr(new Lattice_element_slice(
            elements.front(), 0, slice_end));
    Lattice_element_slices slices;
    slices.push_back(slice_sptr);
    chef_lattice_sptr->construct_sliced_beamline(slices);
    Independent_operations ops = mixed_chef_o_e.extract(b.reference_particle,
            slices);
    BOOST_CHECK_EQUAL(ops.size(), 1);
    BOOST_CHECK_EQUAL(ops.front()->get_type(), fast_mapping_type_name);
}
// test_note: we still need to check the extracted value

BOOST_FIXTURE_TEST_CASE(extract3, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_mixed_operation_extractor mixed_chef_o_e(chef_lattice_sptr, map_order);
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
    Independent_operations ops = mixed_chef_o_e.extract(b.reference_particle,
            slices);
    BOOST_CHECK_EQUAL(ops.size(), 1);
    BOOST_CHECK_EQUAL(ops.front()->get_type(), fast_mapping_type_name);
}
// test_note: we still need to check the extracted value

struct Rf_lattice_fixture
{
    Rf_lattice_fixture() :
        b(), lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup rf lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element r("rfcavity", "r");
        r.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);

        lattice_sptr->append(f);
        lattice_sptr->append(r);
        lattice_sptr->append(d);
        lattice_sptr->append(r);
        lattice_sptr->set_reference_particle(b.reference_particle);
    }
    ;

    ~Rf_lattice_fixture()
    {
        BOOST_TEST_MESSAGE("teardown rf lattice fixture");

    }
    ;

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;
};

BOOST_FIXTURE_TEST_CASE(extract1rf, Rf_lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_mixed_operation_extractor mixed_chef_o_e(chef_lattice_sptr, map_order);
    Lattice_elements elements(lattice_sptr->get_elements());
    Lattice_elements::iterator it = elements.begin();
    ++it;
    Lattice_element_slice_sptr slice_sptr(new Lattice_element_slice(*it));
    Lattice_element_slices slices;
    slices.push_back(slice_sptr);
    chef_lattice_sptr->construct_sliced_beamline(slices);
    Independent_operations ops = mixed_chef_o_e.extract(b.reference_particle,
            slices);
    BOOST_CHECK_EQUAL(ops.size(), 3);
    Independent_operations::const_iterator op_it = ops.begin();
    BOOST_CHECK_EQUAL((*op_it)->get_type(), fast_mapping_type_name);
    ++op_it;
    BOOST_CHECK_EQUAL((*op_it)->get_type(), chef_propagate_type_name);
    ++op_it;
    BOOST_CHECK_EQUAL((*op_it)->get_type(), fast_mapping_type_name);
}
// test_note: we still need to check the extracted value

struct Thinrf_lattice_fixture
{
    Thinrf_lattice_fixture() :
        b(), lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup thinrf lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element o1("drift", "o1");
        o1.set_double_attribute("l", drift_length * 0.5);
        Lattice_element r("rfcavity", "r");
        r.set_double_attribute("l", 0);
        Lattice_element o2("drift", "o2");
        o2.set_double_attribute("l", drift_length * 0.5);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);

        lattice_sptr->append(f);
        lattice_sptr->append(o1);
        lattice_sptr->append(r);
        lattice_sptr->append(o2);
        lattice_sptr->append(d);
        lattice_sptr->append(o1);
        lattice_sptr->append(r);
        lattice_sptr->append(o2);
        lattice_sptr->set_reference_particle(b.reference_particle);
    }
    ;

    ~Thinrf_lattice_fixture()
    {
        BOOST_TEST_MESSAGE("teardown thinrf lattice fixture");

    }
    ;

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;
};

BOOST_FIXTURE_TEST_CASE(extract4thinrf, Thinrf_lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_mixed_operation_extractor mixed_chef_o_e(chef_lattice_sptr, map_order);
    Lattice_elements elements(lattice_sptr->get_elements());

    Lattice_element_slices slices;
    Lattice_elements::iterator it = elements.begin();
    // get four slices...

    // ... the whole first element
    Lattice_element_slice_sptr slice_sptr1(new Lattice_element_slice(*it));
    slices.push_back(slice_sptr1);

    // ... the whole second element
    ++it;
    Lattice_element_slice_sptr slice_sptr2(new Lattice_element_slice(*it));
    slices.push_back(slice_sptr2);

    // ... the whole third element
    ++it;
    Lattice_element_slice_sptr slice_sptr3(new Lattice_element_slice(*it));
    slices.push_back(slice_sptr3);

    /// ... part of the fourth element
    ++it;
    double fourth_fraction = 0.3;
    double fourth_length = (*it)->get_length();
    double fourth_end = fourth_fraction * fourth_length;
    Lattice_element_slice_sptr slice_sptr4(new Lattice_element_slice(*it, 0,
            fourth_end));

    slices.push_back(slice_sptr4);
    chef_lattice_sptr->construct_sliced_beamline(slices);
    Independent_operations ops = mixed_chef_o_e.extract(b.reference_particle,
            slices);
    BOOST_CHECK_EQUAL(ops.size(), 3);
    Independent_operations::iterator op_it = ops.begin();
    BOOST_CHECK_EQUAL((*op_it)->get_type(), fast_mapping_type_name);
    ++op_it;
    BOOST_CHECK_EQUAL((*op_it)->get_type(), chef_propagate_type_name);
    ++op_it;
    BOOST_CHECK_EQUAL((*op_it)->get_type(), fast_mapping_type_name);
}
// test_note: we still need to check the extracted value

BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Chef_lattice_sptr chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    Chef_mixed_operation_extractor mixed_chef_o_e(chef_lattice_sptr, map_order);
    xml_save(mixed_chef_o_e, "chef_mixed_operation_extractor.xml");

    Chef_mixed_operation_extractor loaded;
    xml_load(loaded, "chef_mixed_operation_extractor.xml");
}
