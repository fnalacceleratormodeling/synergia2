#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"



BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-14;

const double mass = 100.0;
const double total_energy = 125.0;
const double real_num = 2.0e12;

struct Fixture
{
    Fixture(int total_num, int total_spectator_num) :
        four_momentum(mass, total_energy), 
        reference_particle(pconstants::proton_charge, four_momentum),
        comm_sptr(new Commxx()), 
        bunch(reference_particle, total_num, total_spectator_num, real_num, comm_sptr),
        total_num(total_num),
        total_spectator_num(total_spectator_num)
    {
        BOOST_TEST_MESSAGE("setup fixture");
    }

    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;

    int total_num;
    int total_spectator_num;
};

struct Fix_p100_s0 : public Fixture
{ Fix_p100_s0() : Fixture(100, 0) { } };

struct Fix_p100_s10 : public Fixture
{ Fix_p100_s10() : Fixture(100, 10) { } };

void
dummy_populate(Bunch &bunch, int offset = 0)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) 
    {
        // coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            bunch.get_local_particles()[part][i] = 10.0 * (part + offset) + i;
        }

        // momenta
        for (int i = 1; i < 6; i += 2) 
        {
            bunch.get_local_particles()[part][i] = 1e-4 * (10.0 * (part + offset) + i);
        }
    }

    // spectator particles
    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        // coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            bunch.get_local_spectator_particles()[part][i] = 10.0 * (part + offset) + i;
        }

        // momenta
        for (int i = 1; i < 6; i += 2) 
        {
            bunch.get_local_spectator_particles()[part][i] = 1e-4 * (10.0 * (part + offset) + i);
        }
    }

}

void
random_populate(Bunch &bunch)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) 
    {
        for (int i = 0; i < 6; i += 1) 
        {
            bunch.get_local_particles()[part][i] = 
                double(std::rand()) / RAND_MAX;
        }

        bunch.get_local_particles()[part][6] = part;
    }

    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        for (int i = 0; i < 6; i += 1) 
        {
            bunch.get_local_spectator_particles()[part][i] = 
                double(std::rand()) / RAND_MAX;
        }

        bunch.get_local_spectator_particles()[part][6] = part;
    }
}

void
compare_bunches(Bunch &bunch1, Bunch &bunch2, double tolerance,
        bool check_state = true, bool check_ids = true, bool check_bucket_index = true)
{
    BOOST_CHECK_EQUAL(
            bunch1.get_reference_particle().get_total_energy(),
            bunch2.get_reference_particle().get_total_energy() );

    BOOST_CHECK_EQUAL(
            bunch1.get_particle_charge(),
            bunch2.get_particle_charge() );

    BOOST_CHECK_CLOSE(bunch1.get_mass(), bunch2.get_mass(), tolerance);
    BOOST_CHECK_CLOSE(bunch1.get_real_num(), bunch2.get_real_num(), tolerance);
    BOOST_CHECK_EQUAL(bunch1.get_local_num(), bunch2.get_local_num());
    BOOST_CHECK_EQUAL(bunch1.get_total_num(), bunch2.get_total_num());

    if (check_bucket_index) 
    {
        BOOST_CHECK_EQUAL(bunch1.is_bucket_index_assigned(), bunch2.is_bucket_index_assigned());
        if (bunch1.is_bucket_index_assigned()){
            BOOST_CHECK_EQUAL(bunch1.get_bucket_index(), bunch2.get_bucket_index());
        }
    }

    BOOST_CHECK_EQUAL(bunch1.get_z_period_length(), bunch2.get_z_period_length());
    BOOST_CHECK_EQUAL(bunch1.get_longitudinal_aperture_length(), bunch2.get_longitudinal_aperture_length());
    BOOST_CHECK_EQUAL(bunch1.is_z_periodic(), bunch2.is_z_periodic());
    BOOST_CHECK_EQUAL(bunch1.has_longitudinal_aperture(),  bunch2.has_longitudinal_aperture());

    if (check_state) 
    {
        BOOST_CHECK_EQUAL(bunch1.get_state(), bunch2.get_state());
    }

    for (int part = 0; part < bunch1.get_local_num(); ++part) 
    {
        // this loop is unrolled in order to give more meaningful error messages
        // i.e., error messages including which component was tested
        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][0],
                bunch2.get_local_particles()[part][0], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][1],
                bunch2.get_local_particles()[part][1], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][2],
                bunch2.get_local_particles()[part][2], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][3],
                bunch2.get_local_particles()[part][3], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][4],
                bunch2.get_local_particles()[part][4], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][5],
                bunch2.get_local_particles()[part][5], tolerance );

        if (check_ids) 
        {
            BOOST_CHECK_CLOSE(
                    bunch1.get_local_particles()[part][6],
                    bunch2.get_local_particles()[part][6], tolerance );
        }
    }

    for (int part = 0; part < bunch1.get_local_spectator_num(); ++part) 
    {
        // this loop is unrolled in order to give more meaningful error messages
        // i.e., error messages including which component was tested
        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][0],
                bunch2.get_local_spectator_particles()[part][0], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][1],
                bunch2.get_local_spectator_particles()[part][1], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][2],
                bunch2.get_local_spectator_particles()[part][2], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][3],
                bunch2.get_local_spectator_particles()[part][3], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][4],
                bunch2.get_local_spectator_particles()[part][4], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][5],
                bunch2.get_local_spectator_particles()[part][5], tolerance );

        if (check_ids) 
        {
            BOOST_CHECK_CLOSE(
                    bunch1.get_local_spectator_particles()[part][6],
                    bunch2.get_local_spectator_particles()[part][6], tolerance );
        }
    }
}


BOOST_FIXTURE_TEST_CASE(construct_no_spectator, Fix_p100_s0)
{
}

BOOST_FIXTURE_TEST_CASE(construct_spectator, Fix_p100_s10)
{
}

BOOST_FIXTURE_TEST_CASE(spectator_num, Fix_p100_s10)
{
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), 10);
    BOOST_CHECK_EQUAL(bunch.get_total_spectator_num(), 10);
}

BOOST_FIXTURE_TEST_CASE(get_local_spectator_particles, Fix_p100_s10)
{
    MArray2d_ref particles(bunch.get_local_spectator_particles());
    BOOST_CHECK_EQUAL(particles.shape()[1], 7);

    unsigned int u_local_num = bunch.get_local_spectator_num();
    BOOST_CHECK(particles.shape()[0] >= u_local_num);

    unsigned int u_local_num_slots = bunch.get_local_spectator_num_slots();
    BOOST_CHECK(particles.shape()[0] == u_local_num_slots);
}

BOOST_FIXTURE_TEST_CASE(get_const_local_spectator_particles, Fix_p100_s10)
{
    Const_MArray2d_ref particles(bunch.get_local_spectator_particles());
    BOOST_CHECK_EQUAL(particles.shape()[1], 7);

    unsigned int u_local_num = bunch.get_local_spectator_num();
    BOOST_CHECK(particles.shape()[0] >= u_local_num);

    unsigned int u_local_num_slots = bunch.get_local_spectator_num_slots();
    BOOST_CHECK(particles.shape()[0] == u_local_num_slots);
}

BOOST_FIXTURE_TEST_CASE(set_local_spectator_num, Fix_p100_s10)
{
    int new_local_num = total_spectator_num / comm_sptr->get_size() - 5;

    bunch.set_local_spectator_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), 100);
}

BOOST_FIXTURE_TEST_CASE(set_local_spectator_num_neg, Fix_p100_s10)
{
    int new_local_num = total_spectator_num / comm_sptr->get_size() - 12;

    bunch.set_local_spectator_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), 0);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), 100);
}

BOOST_FIXTURE_TEST_CASE(get_total_spectator_num, Fix_p100_s10)
{
    BOOST_CHECK_EQUAL(bunch.get_total_spectator_num(), 10);
}

BOOST_FIXTURE_TEST_CASE(update_total_num_with_new_spectator, Fix_p100_s10)
{
    const int new_local_num = bunch.get_local_num() - 7;
    bunch.set_local_spectator_num(new_local_num);
    bunch.update_total_num();

    BOOST_CHECK_EQUAL(
            bunch.get_total_spectator_num(), 
            new_local_num*Commxx().get_size() );

    BOOST_CHECK_EQUAL(bunch.get_total_num(), 100);
}

#if 0
BOOST_FIXTURE_TEST_CASE(set_local_num, Fix_p100_s10)
{
    int new_local_num = total_num / comm_sptr->get_size() - 5;

    bunch.set_local_spectator_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), 100);
}

BOOST_FIXTURE_TEST_CASE(get_total_spectator_num, Fix_p100_s10)
{
    BOOST_CHECK_EQUAL(bunch.get_total_spectator_num(), 10);
}

BOOST_FIXTURE_TEST_CASE(update_total_num_with_new_spectator, Fix_p100_s10)
{
    const int new_local_num = bunch.get_local_num() - 7;
    bunch.set_local_spectator_num(new_local_num);
    bunch.update_total_num();

    BOOST_CHECK_EQUAL(
            bunch.get_total_spectator_num(), 
            new_local_num*Commxx().get_size() );

    BOOST_CHECK_EQUAL(bunch.get_total_num(), 100);
}
#endif



#if 0
BOOST_FIXTURE_TEST_CASE(set_bucket_index, Fixture)
{
  bunch.set_bucket_index(4);
  BOOST_CHECK_EQUAL(bunch.is_bucket_index_assigned(),true);
  BOOST_CHECK_EQUAL(bunch.get_bucket_index(),4);
}

BOOST_FIXTURE_TEST_CASE(set_z_period, Fixture)
{
  double z_length=5.5;
  bunch.set_z_period_length(z_length);
  BOOST_CHECK_EQUAL(bunch.is_z_periodic(),true);
  BOOST_CHECK_EQUAL(bunch.has_longitudinal_aperture(), false);
  BOOST_CHECK_EQUAL(bunch.get_z_period_length(),z_length);
}

BOOST_FIXTURE_TEST_CASE(set_longitudinal_aperture, Fixture)
{
  double z_length=5.5;
  bunch.set_longitudinal_aperture_length(z_length);
  BOOST_CHECK_EQUAL(bunch.is_z_periodic(),false);
  BOOST_CHECK_EQUAL(bunch.has_longitudinal_aperture(), true);
  BOOST_CHECK_EQUAL(bunch.get_longitudinal_aperture_length(),z_length);
}

BOOST_FIXTURE_TEST_CASE(conflict_longitudinal_aperture_z_periodic, Fixture)
{
  double z_length=5.5;
  bunch.set_longitudinal_aperture_length(z_length);
  bool caught_error = false;
  try {
      bunch.set_z_period_length(z_length);
  }catch (std::runtime_error) {
        caught_error = true;
  }
  BOOST_CHECK(caught_error);
  
}

BOOST_FIXTURE_TEST_CASE(conflict_z_periodic_longitudinal_aperture, Fixture)
{
  double z_length=5.5;
  bunch.set_z_period_length(z_length);
  bool caught_error = false;
  try {      
     bunch.set_longitudinal_aperture_length(z_length);
  }catch (std::runtime_error) {
        caught_error = true;
  }
  BOOST_CHECK(caught_error);
  
} 

BOOST_FIXTURE_TEST_CASE(check_ids, Fixture)
{
    double offset = bunch.get_local_particles()[0][Bunch::id];
    if (comm_sptr->get_rank() == 0) {
        for (int part = 0; part < bunch.get_local_num(); ++part) {
            BOOST_CHECK_CLOSE(part + offset, bunch.get_local_particles()[part][Bunch::id], tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(check_ids2, Fixture)
{
    Bunch second_bunch(reference_particle, total_num, real_num, comm_sptr);
    BOOST_CHECK( static_cast<int>(second_bunch.get_local_particles()[0][Bunch::id]
                    - bunch.get_local_particles()[0][Bunch::id]) == total_num);
}

BOOST_FIXTURE_TEST_CASE(copy_construct, Fixture)
{
    dummy_populate(bunch);
    Bunch second_bunch(bunch);
    compare_bunches(bunch, second_bunch, tolerance);
}

BOOST_FIXTURE_TEST_CASE(assign, Fixture)
{
    Bunch second_bunch(reference_particle, total_num + 10, real_num * 2, comm_sptr);
    dummy_populate(bunch);
    second_bunch = bunch;
    compare_bunches(bunch, second_bunch, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_particle_charge, Fixture)
{
    BOOST_CHECK_EQUAL(bunch.get_particle_charge(), pconstants::proton_charge);
}


BOOST_FIXTURE_TEST_CASE(get_mass, Fixture)
{
    BOOST_CHECK_CLOSE(bunch.get_mass(), mass, tolerance);
}

BOOST_FIXTURE_TEST_CASE(set_particle_charge, Fixture)
{
    bunch.set_particle_charge(pconstants::electron_charge);
    BOOST_CHECK_EQUAL(bunch.get_particle_charge(), pconstants::electron_charge);
}

BOOST_FIXTURE_TEST_CASE(get_real_num, Fixture)
{
    BOOST_CHECK_CLOSE(bunch.get_real_num(), real_num, tolerance);
}

BOOST_FIXTURE_TEST_CASE(set_real_num, Fixture)
{
    const double new_real_num = real_num * 1.5;
    bunch.set_real_num(new_real_num);
    BOOST_CHECK_CLOSE(bunch.get_real_num(), new_real_num, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_local_num, Fixture)
{
    // n.b.: this test assumes that we are running on one processor
    BOOST_CHECK_EQUAL(bunch.get_local_num(), total_num);
}

BOOST_FIXTURE_TEST_CASE(set_local_num, Fixture)
{
    int new_local_num = total_num / comm_sptr->get_size() - 5;
    bunch.set_local_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), new_local_num);
}

BOOST_FIXTURE_TEST_CASE(get_total_num, Fixture)
{
    BOOST_CHECK_EQUAL(bunch.get_total_num(), total_num);
}

BOOST_FIXTURE_TEST_CASE(update_total_num, Fixture)
{
    const int new_local_num = bunch.get_local_num() - 7;
    bunch.set_local_num(new_local_num);
    bunch.update_total_num();
    BOOST_CHECK_EQUAL(bunch.get_total_num(),
            new_local_num*Commxx().get_size());
}

BOOST_FIXTURE_TEST_CASE(set_get_sort_period, Fixture)
{
    const int new_period = 17;
    bunch.set_sort_period(new_period);
    BOOST_CHECK_EQUAL(bunch.get_sort_period(), new_period);
}

bool
sorted(Bunch const& bunch, int index)
{
    bool retval = true;
    for (int part = 1; part < bunch.get_local_num(); ++part) {
        if (bunch.get_local_particles()[part][index]
                < bunch.get_local_particles()[part - 1][index]) {
            retval = false;
        }
    }
    return retval;
}

BOOST_FIXTURE_TEST_CASE(sort, Fixture)
{
    random_populate(bunch);
    for (int sort_index = 0; sort_index < 6; ++sort_index) {
        bunch.sort(sort_index);
        BOOST_CHECK(sorted(bunch, sort_index));
    }
}

BOOST_FIXTURE_TEST_CASE(sort_bad, Fixture)
{
    bool caught_error = false;
    try {
        bunch.sort(-1);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);

    caught_error = false;
    try {
        bunch.sort(7);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(periodic_sort, Fixture)
{
    random_populate(bunch);
    const int sort_index = 1;
    const int sort_period = 7;
    bunch.set_sort_period(sort_period);
    for (int i = 0; i < sort_period; ++i) {
        bunch.periodic_sort(sort_index);
        BOOST_CHECK(!sorted(bunch, sort_index));
    }
    bunch.periodic_sort(sort_index);
    BOOST_CHECK(sorted(bunch, sort_index));
}

BOOST_FIXTURE_TEST_CASE(get_reference_particle, Fixture)
{
    Reference_particle ref(bunch.get_reference_particle());
}

BOOST_FIXTURE_TEST_CASE(get_const_reference_particle, Fixture)
{
    const Reference_particle ref(bunch.get_reference_particle());
}

BOOST_FIXTURE_TEST_CASE(get_local_particles, Fixture)
{
    MArray2d_ref local_particles(bunch.get_local_particles());
    BOOST_CHECK_EQUAL(local_particles.shape()[1], 7);
    unsigned int u_local_num = bunch.get_local_num();
    BOOST_CHECK(local_particles.shape()[0] >= u_local_num);
}

BOOST_FIXTURE_TEST_CASE(get_const_local_particles, Fixture)
{
    Const_MArray2d_ref local_particles(bunch.get_local_particles());
    BOOST_CHECK_EQUAL(local_particles.shape()[1], 7);
    unsigned int u_local_num = bunch.get_local_num();
    BOOST_CHECK(local_particles.shape()[0] >= u_local_num);
}

BOOST_FIXTURE_TEST_CASE(increase_local_num, Fixture)
{
    const int increase = 5;
    int old_num = bunch.get_local_num();
    bunch.set_local_num(old_num + increase);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), old_num + increase);
    // make sure I can scribble in the new space
    for (int part=old_num; part<bunch.get_local_num(); ++part) {
        for (int i=0; i<6; ++i) {
            bunch.get_local_particles()[part][i] = -100.0;
        }
        bunch.get_local_particles()[part][6] = part;
    }
}

BOOST_FIXTURE_TEST_CASE(get_state, Fixture)
{
    Bunch::State state;
    state = bunch.get_state();
    BOOST_CHECK_EQUAL(state,Bunch::fixed_z_lab);
}

BOOST_FIXTURE_TEST_CASE(get_comm, Fixture)
{
    BOOST_CHECK_EQUAL(bunch.get_comm().get(),MPI_COMM_WORLD);
}

BOOST_FIXTURE_TEST_CASE(convert_to_state, Fixture)
{
    dummy_populate(bunch);
    // shift the bunch to have zero mean in all coordinates
    MArray1d orig_bunch_mean(Core_diagnostics::calculate_mean(bunch));
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            bunch.get_local_particles()[part][i] -= orig_bunch_mean[i];
        }
    }

    Bunch second_bunch(bunch);
    bunch.convert_to_state(Bunch::fixed_t_bunch);

    BOOST_CHECK(bunch.get_state() == Bunch::fixed_t_bunch);
    BOOST_CHECK(second_bunch.get_state() == Bunch::fixed_z_lab);

    // verify that transverse means are unaffected by state transformation
    MArray1d accel_mean(Core_diagnostics::calculate_mean(second_bunch));
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    const double tolerance_mean = 1.0e-11;
    BOOST_CHECK(floating_point_equal(accel_mean[Bunch::x], bunch_mean[Bunch::x], tolerance_mean));
    BOOST_CHECK(floating_point_equal(accel_mean[Bunch::xp], bunch_mean[Bunch::xp], tolerance_mean));
    BOOST_CHECK(floating_point_equal(accel_mean[Bunch::y], bunch_mean[Bunch::y], tolerance_mean));
    BOOST_CHECK(floating_point_equal(accel_mean[Bunch::yp], bunch_mean[Bunch::yp], tolerance_mean));

    // verify that transverse std's are unaffected by state transformation
    MArray1d accel_std(Core_diagnostics::calculate_std(second_bunch, accel_mean));
    MArray1d bunch_std(Core_diagnostics::calculate_std(bunch, bunch_mean));
    const double tolerance_std = 1.0e-11;
    BOOST_CHECK_CLOSE(accel_std[Bunch::x], bunch_std[Bunch::x], tolerance_std);
    BOOST_CHECK_CLOSE(accel_std[Bunch::xp], bunch_std[Bunch::xp], tolerance_std);
    BOOST_CHECK_CLOSE(accel_std[Bunch::y], bunch_std[Bunch::y], tolerance_std);
    BOOST_CHECK_CLOSE(accel_std[Bunch::yp], bunch_std[Bunch::yp], tolerance_std);

    // verify that longitudinal std's are transformed as expected
    double beta = bunch.get_reference_particle().get_beta();
    double gamma = bunch.get_reference_particle().get_gamma();
    BOOST_CHECK_CLOSE(gamma * beta * accel_std[Bunch::cdt],
            bunch_std[Bunch::z], tolerance_std);

    // Check relativistic invariants for spatial vectors
    const double tolerance_invariant = 1.0e-15;
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x, y, z, cdt;
        // accelerator frame
        x = second_bunch.get_local_particles()[part][Bunch::x];
        y = second_bunch.get_local_particles()[part][Bunch::y];
        z = 0;
        cdt = second_bunch.get_local_particles()[part][Bunch::cdt];
        double accel_invariant = x * x + y * y + z * z - cdt * cdt;
        // bunch frame
        x = bunch.get_local_particles()[part][Bunch::x];
        y = bunch.get_local_particles()[part][Bunch::y];
        z = bunch.get_local_particles()[part][Bunch::z];
        cdt = gamma * second_bunch.get_local_particles()[part][Bunch::cdt];
        double bunch_invariant = x * x + y * y + z * z - cdt * cdt;
        BOOST_CHECK_CLOSE(accel_invariant, bunch_invariant, tolerance_invariant);
    }

    // Check that the roundtrip transformation is the identity
    bunch.convert_to_state(Bunch::fixed_z_lab);
    const double convert_tolerance = 5.0e-8;
    compare_bunches(bunch, second_bunch, convert_tolerance);
}

BOOST_FIXTURE_TEST_CASE(convert_to_fixed_t_bad_pz, Fixture)
{
    dummy_populate(bunch);
    const int bad_particle_id = total_num / 2;
    bunch.get_local_particles()[bad_particle_id][Bunch::xp] = 0.8;
    bunch.get_local_particles()[bad_particle_id][Bunch::yp] = -0.7;
    bool caught_error = false;
    try {
        bunch.convert_to_state(Bunch::fixed_t_bunch);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

class Fixed_t_z_dummy : public Fixed_t_z_converter
{
public:
 //   void
 //   fixed_t_to_fixed_z(Bunch &bunch)
 //   {
 //   }
 //   ;
 //   void
 //   fixed_z_to_fixed_t(Bunch &bunch)
 //   {
//    }
  //  ;



    void
    from_z_lab_to_t_lab(Bunch &bunch){};

    void
    from_t_lab_to_z_lab(Bunch &bunch){};

    void
    from_z_lab_to_t_bunch(Bunch &bunch){};

    void
    from_t_bunch_to_z_lab(Bunch &bunch){};

    void
    from_t_lab_to_t_bunch(Bunch &bunch){};

    void
    from_t_bunch_to_t_lab(Bunch &bunch){};

};

BOOST_FIXTURE_TEST_CASE(set_converter, Fixture)
{
    // This test relies on the Fixed_t_z_dummy class not
    // doing anything to the bunch. It verifies that we are not
    // using the default converter after set_converter.
    Fixed_t_z_dummy converter;
    bunch.set_converter(converter);
    dummy_populate(bunch);
    Bunch second_bunch(bunch);
    bunch.convert_to_state(Bunch::fixed_t_bunch);
    compare_bunches(bunch, second_bunch, tolerance, false);
}

// inject tolerance is lower than the others because calculation of
// momentum rescaling induces roundoff errors
const double inject_tolerance = 2.0e11;

BOOST_FIXTURE_TEST_CASE(inject, Fixture)
{
    Bunch total_bunch(bunch);
    const int local_num = 10;
    bunch.set_local_num(local_num);
    bunch.update_total_num();
    Bunch second_bunch(bunch);
    dummy_populate(bunch);
    dummy_populate(second_bunch, local_num);
    total_bunch.set_local_num(2 * local_num);
    total_bunch.update_total_num();
    dummy_populate(total_bunch);
    bunch.inject(second_bunch);
    compare_bunches(bunch, total_bunch, inject_tolerance, true, false);
}

BOOST_FIXTURE_TEST_CASE(inject2, Fixture)
{
    Bunch total_bunch(bunch);
    const int local_num = bunch.get_local_num();
    Bunch second_bunch(bunch);
    dummy_populate(second_bunch);
    dummy_populate(total_bunch);
    total_bunch.inject(second_bunch);
    // total bunch should now have two copies of bunch
    BOOST_CHECK_EQUAL(total_bunch.get_local_num(), 2*local_num);
    BOOST_CHECK_EQUAL(total_bunch.get_total_num(), 2*bunch.get_total_num());
    BOOST_CHECK_CLOSE(total_bunch.get_real_num(), 2.0*bunch.get_real_num(), tolerance);
    for (int part=0; part<local_num; ++part) {
        for (int coord=0; coord<6; ++coord) {
            BOOST_CHECK_CLOSE(total_bunch.get_local_particles()[part][coord],
                              total_bunch.get_local_particles()[part+local_num][coord], inject_tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(inject_mismatched_weights, Fixture)
{
    Bunch second_bunch(bunch);
    second_bunch.set_real_num(second_bunch.get_real_num() * 2.0);
    bool caught_error = false;
    try {
        bunch.inject(second_bunch);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(inject_wrong_particle, Fixture)
{
    Four_momentum four_momentum2(mass/2.0, total_energy);
    Reference_particle reference_particle2(1, four_momentum2);
    Bunch target_bunch(bunch);
    Bunch second_bunch(reference_particle2, total_num, real_num, comm_sptr);
    bool caught_error = false;

    try {
        target_bunch.inject(second_bunch);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);

    Bunch target_bunch2(bunch);
    Reference_particle reference_particle3(2, reference_particle.get_four_momentum());
    Bunch third_bunch(reference_particle3, total_num, real_num, comm_sptr);
    caught_error = false;
    try {
        target_bunch2.inject(third_bunch);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);

}

BOOST_AUTO_TEST_CASE(inject_different_momentum_bunch)
{
    const double momentum1 = 10.0;
    const double momentum2 = 10.1;
    const double mass = 0.75;
    Four_momentum four_momentum1(mass);
    Four_momentum four_momentum2(mass);
    four_momentum1.set_momentum(momentum1);
    Reference_particle refpart1(1, four_momentum1);
    four_momentum2.set_momentum(momentum2);
    Reference_particle refpart2(1, four_momentum2);
    Commxx_sptr comm_sptr(new Commxx());

    // first bunch has no particles to make it easy
    Bunch bunch1(refpart1, 0, 0.0, comm_sptr);
    Bunch bunch2(refpart2, 3, 1.0e4, comm_sptr);
    // give some transverse momentum to the bunch 2 particles
    for (int part=0; part<bunch2.get_local_num(); ++part) {
        bunch2.get_local_particles()[part][1] = .001*(part-1);
        bunch2.get_local_particles()[part][3] = -.001*(part-1);
        // try different longitudinal momenta
        bunch2.get_local_particles()[part][5] = 0.01*(part-1);
    }
    // inject it
    bunch1.inject(bunch2);
    // the unscaled momentum better be the same
    for (int part=0; part<bunch1.get_local_num(); ++part) {
        BOOST_CHECK_CLOSE(momentum1*bunch1.get_local_particles()[part][1],
                          momentum2*bunch2.get_local_particles()[part][1], inject_tolerance);
        BOOST_CHECK_CLOSE(momentum1*bunch1.get_local_particles()[part][3],
                          momentum2*bunch2.get_local_particles()[part][3], inject_tolerance);
        BOOST_CHECK_CLOSE(momentum1*(1.0+bunch1.get_local_particles()[part][5]),
                          momentum2*(1.0+bunch2.get_local_particles()[part][5]), inject_tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(inject_into_empty_bunch, Fixture)
{
    Bunch new_bunch(reference_particle, 0, 0.0, comm_sptr);
    Bunch bunch_to_inject(bunch);
    dummy_populate(bunch_to_inject);
    new_bunch.inject(bunch_to_inject);
    // new_bunch should be the same as old bunch
    compare_bunches(new_bunch, bunch_to_inject, inject_tolerance, true, false, false);
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Fixture)
{
    dummy_populate(bunch);
    xml_save(bunch, "bunch.xml");

    Bunch loaded;
    xml_load(loaded, "bunch.xml");

    compare_bunches(bunch, loaded, tolerance);
}

BOOST_FIXTURE_TEST_CASE(default_construct_no_segfault, Fixture)
{
    Bunch b;
}
#endif

