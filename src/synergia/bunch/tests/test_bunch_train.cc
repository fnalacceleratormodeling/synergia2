#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/bunch_train.h"
#include "bunches_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double default_tolerance = 1.0e-14;

void
compare_bunches(Bunch &bunch1, Bunch &bunch2, double tolerance = default_tolerance,
        bool check_state = true, bool check_ids = true)
{
    BOOST_CHECK_EQUAL(bunch1.get_reference_particle().get_total_energy(),
            bunch2.get_reference_particle().get_total_energy());
    BOOST_CHECK_EQUAL(bunch1.get_particle_charge(),
            bunch2.get_particle_charge());
    BOOST_CHECK_CLOSE(bunch1.get_mass(), bunch2.get_mass(), tolerance);
    BOOST_CHECK_CLOSE(bunch1.get_real_num(), bunch2.get_real_num(), tolerance);
    BOOST_CHECK_EQUAL(bunch1.get_local_num(), bunch2.get_local_num());
    BOOST_CHECK_EQUAL(bunch1.get_total_num(), bunch2.get_total_num());
    BOOST_CHECK_EQUAL(bunch1.is_bucket_index_assigned(), bunch2.is_bucket_index_assigned());
    if (bunch1.is_bucket_index_assigned()){
        BOOST_CHECK_EQUAL(bunch1.get_bucket_index(), bunch2.get_bucket_index());
    }

    BOOST_CHECK_EQUAL(bunch1.get_longitudinal_boundary().first, bunch2.get_longitudinal_boundary().first);
    BOOST_CHECK_EQUAL(bunch1.get_longitudinal_boundary().second, bunch2.get_longitudinal_boundary().second);


    if (check_state) {
        BOOST_CHECK_EQUAL(bunch1.get_state(), bunch2.get_state());
    }
    for (int part = 0; part < bunch1.get_local_num(); ++part) {
        // this loop is unrolled in order to give more meaningful error messages
        // i.e., error messages including which component was tested
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][0],
                bunch2.get_local_particles()[part][0], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][1],
                bunch2.get_local_particles()[part][1], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][2],
                bunch2.get_local_particles()[part][2], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][3],
                bunch2.get_local_particles()[part][3], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][4],
                bunch2.get_local_particles()[part][4], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][5],
                bunch2.get_local_particles()[part][5], tolerance);
        if (check_ids) {
            BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][6],
                    bunch2.get_local_particles()[part][6], tolerance);
        }
    }
}



BOOST_FIXTURE_TEST_CASE(construct1, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);
}

BOOST_FIXTURE_TEST_CASE(construct2, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);
}

BOOST_FIXTURE_TEST_CASE(construct3, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    const double bunch_separation3 = 99.9;
    separations.push_back(bunch_separation3);
    bool caught_error = false;
    try {
        Bunch_train bunch_train(bunches, separations);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(get_num_bunches, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);

    BOOST_CHECK_EQUAL(bunch_train.get_size(), num_bunches);
}

BOOST_FIXTURE_TEST_CASE(get_bunches, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);

    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(0), bunches.at(0));
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(1), bunches.at(1));
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(2), bunches.at(2));
}

BOOST_FIXTURE_TEST_CASE(set_bucket_indices1, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);

    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(0)->get_bucket_index(), 0);
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(1)->get_bucket_index(), 1);
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(2)->get_bucket_index(), 2);
}

BOOST_FIXTURE_TEST_CASE(set_bucket_indices2, Bunches_fixture)
{
    bunches.at(0)->set_bucket_index(1);
    bunches.at(1)->set_bucket_index(7);
    bunches.at(2)->set_bucket_index(83);
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);

    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(0)->get_bucket_index(), 1);
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(1)->get_bucket_index(), 7);
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(2)->get_bucket_index(), 83);
}

BOOST_FIXTURE_TEST_CASE(set_bucket_indices_bad, Bunches_fixture)
{
    bunches.at(0)->set_bucket_index(0);
    bunches.at(1)->set_bucket_index(83);
    bunches.at(2)->set_bucket_index(83);
    const double bunch_separation = 1.7;
    bool caught = false;
    try {
        Bunch_train bunch_train(bunches, bunch_separation);
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_FIXTURE_TEST_CASE(get_spacings1, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);

    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(0), bunch_separation1);
    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(1), bunch_separation2);
}

BOOST_FIXTURE_TEST_CASE(get_spacings2, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);

    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(0), bunch_separation);
    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(1), bunch_separation);
}

BOOST_FIXTURE_TEST_CASE(get_parent_comm_sptr, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);    
    int result;
    MPI_Comm_compare(MPI_COMM_WORLD, bunch_train.get_parent_comm_sptr()->get(), &result);
    BOOST_CHECK(result == MPI_IDENT);    
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);
    xml_save(bunch_train, "bunch_train.xml");
    
    Bunch_train bunch_loaded;
    xml_load(bunch_loaded, "bunch_train.xml");
    size_t num_bunches=bunch_train.get_bunches().size();
    size_t num_loaded=bunch_loaded.get_bunches().size();
    BOOST_CHECK_EQUAL(num_bunches, num_loaded);

    int result;
    MPI_Comm_compare( bunch_loaded.get_parent_comm_sptr()->get(),
		      bunch_train.get_parent_comm_sptr()->get(), &result);
    BOOST_CHECK(result == MPI_IDENT);    	      
    
    for (size_t i=0; i<num_bunches; ++i){
      compare_bunches( *bunch_train.get_bunches().at(i), *bunch_loaded.get_bunches().at(i));
    }
    
 
}
