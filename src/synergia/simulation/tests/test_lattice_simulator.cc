#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/lattice_simulator.h"
#include "lattice_fixture.h"
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
}

BOOST_FIXTURE_TEST_CASE(set_slices, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*it, 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*it, 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.set_slices(slices);
}

BOOST_FIXTURE_TEST_CASE(get_map_order, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BOOST_CHECK_EQUAL(lattice_simulator.get_map_order(), map_order);
}

BOOST_FIXTURE_TEST_CASE(get_operation_extractor_map_sptr, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    std::list < std::string
            > names(
                    lattice_simulator.get_operation_extractor_map_sptr()->get_extractor_names());
    std::list < std::string > expected_names;
    expected_names.push_back(default_operation_extractor_name);
    expected_names.push_back(chef_mixed_operation_extractor_name);
    expected_names.push_back(chef_map_operation_extractor_name);
    expected_names.push_back(chef_propagate_operation_extractor_name);

    BOOST_CHECK_EQUAL(names.size(), expected_names.size());
    names.sort();
    expected_names.sort();
    for (std::list<std::string >::iterator it = names.begin(), expected_it =
            expected_names.begin(); it != names.end(); ++it, ++expected_it) {
        BOOST_CHECK_EQUAL((*it), (*expected_it));
    }

    BOOST_CHECK_EQUAL(lattice_simulator.get_map_order(), map_order);
}

BOOST_FIXTURE_TEST_CASE(get_lattice, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.get_lattice();
}

BOOST_FIXTURE_TEST_CASE(get_lattice_sptr, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.get_lattice_sptr();
}

BOOST_FIXTURE_TEST_CASE(get_chef_lattice, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.get_chef_lattice();
}

BOOST_FIXTURE_TEST_CASE(get_chef_lattice_sptr, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.get_chef_lattice_sptr();
}

BOOST_AUTO_TEST_CASE(update)
{
    const double quad_length = 0.2;
    const double quad_strength = 3.2;
    const double drift_length = 3.0;

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    f.set_double_attribute("k1", quad_strength);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", quad_strength);

    Lattice_sptr lattice_sptr(new Lattice(name));
    lattice_sptr->append(f);
    lattice_sptr->append(o);
    lattice_sptr->append(d);
    lattice_sptr->append(o);

    const int charge = pconstants::proton_charge;
    const double mass = pconstants::mp;
    const double total_energy = 125.0;
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*it, 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*it, 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.set_slices(slices);

    double orig_quad_strength = 0.0;
    for (beamline::deep_iterator it =
            lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_begin();
            it
                    != lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_end();
            ++it) {
        if (std::string((*it)->Type()) == "quadrupole") {
            orig_quad_strength = (*it)->Strength();
        }
    }

    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "quadrupole") {
            (*it)->set_double_attribute("k1", 2 * quad_strength);
        }
    }

    lattice_simulator.update();

    double new_quad_strength;
    for (beamline::deep_iterator it =
            lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_begin();
            it
                    != lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_end();
            ++it) {
        if (std::string((*it)->Type()) == "quadrupole") {
            new_quad_strength = (*it)->Strength();
        }
    }
    BOOST_CHECK_CLOSE(new_quad_strength, 2*orig_quad_strength, tolerance);
}

#define CO_DEBUG 0
BOOST_AUTO_TEST_CASE(get_closed_orbit)
{
    const double bend_length = 2.0;
    const double quad_length = 0.5;
    const double sep = 10.0;
    const double focus = 7.0;
    const double quad_strength = 1.0/(focus*quad_length);
    const double drift_length = (sep - quad_length - bend_length)/2.0;
    const double kick_strength = 0.0002;
    const int ncells = 128;
    // angle = 2*pi/(2*ncells)
    const double bend_angle = mconstants::pi/ncells;
    const double co_tolerance = 1.0e-4;

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    f.set_double_attribute("k1", quad_strength);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", -quad_strength);
    Lattice_element b("sbend", "b");
    b.set_double_attribute("l", bend_length);
    b.set_double_attribute("angle", bend_angle);

    Lattice_element k("hkicker", "k");
    k.set_double_attribute("kick", kick_strength);

    Lattice_sptr lattice_sptr(new Lattice(name));
    for (int cell=0; cell<ncells; ++cell) {
        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
    }
    lattice_sptr->append(k);

    const int charge = pconstants::proton_charge;
    const double mass = pconstants::mp;
    const double total_energy = 2.0;
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BmlPtr beamline_sptr(lattice_simulator.get_chef_lattice_sptr()->get_beamline_sptr());
    Proton probe(total_energy);

#if CO_DEBUG
    beamline_sptr->propagate(probe);

    std::cout << "egs: test_get_closed_orbit propagate zero particle: ";
    std::cout << probe.get_x() << " " << probe.get_npx() << " " << probe.get_y() << " " << probe.get_npy()
              << " " << probe.get_cdt() << " " << probe.get_ndp();
    std::cout << std::endl;
#endif

    MArray1d closed_orbit(lattice_simulator.get_closed_orbit());
#if CO_DEBUG
    std::cout << "egs: test_get_closed_orbit: ";
    for (int i=0; i<6; ++i) {
        std::cout << " " << closed_orbit[i];
    }
    std::cout << std::endl;
#endif

    // with kick, the closed orbit in x and xp better be away from 0
    for (int i=0; i<2; ++i) {
        BOOST_CHECK(std::abs(closed_orbit[i]) > 1.0e-5);
    }

    probe.set_x(closed_orbit[0]);
    probe.set_npx(closed_orbit[1]);
    probe.set_y(closed_orbit[2]);
    probe.set_npy(closed_orbit[3]);
    probe.set_cdt(closed_orbit[4]);
    probe.set_ndp(closed_orbit[5]);

    beamline_sptr->propagate(probe);

#if CO_DEBUG
    std::cout << "egs: test_get_closed_orbit after propagate: ";
    std::cout << probe.get_x() << " " << probe.get_npx() << " " << probe.get_y() << " " << probe.get_npy()
              << " " << probe.get_cdt() << " " << probe.get_ndp();
    std::cout << std::endl;
#endif

    BOOST_CHECK(floating_point_equal(closed_orbit[0], probe.get_x(), co_tolerance));
    BOOST_CHECK(floating_point_equal(closed_orbit[1], probe.get_npx(), co_tolerance));
    BOOST_CHECK(floating_point_equal(closed_orbit[2], probe.get_y(), co_tolerance));
    BOOST_CHECK(floating_point_equal(closed_orbit[3], probe.get_npy(), co_tolerance));
}


// void map_applied(MArray2d const& map, MArray1d const& clo)
// {
//   
//     MArray1d  retval(boost::extents[6]);
//      for (int i = 0; i < 6; ++i) {
//        retval[i]=0.;
//        for (int j = 0; j < 6; ++j) {
//           retval[i] += map[i][j]*clo[j];
//        }
//      }
//    
//     
//     std::cout << "coordinates x,xp,y,yp,cdt,dpop, after mapping are[6] \n";
//      for (int i = 0; i < 6; ++i) {
//         std::cout << "x["<<i<<"]="<<" "<<retval[i]<<std::endl;
//      }
//      std::cout << "\n";    
// }

BOOST_AUTO_TEST_CASE(register_closed_orbit)
{
    const double bend_length = 2.0;
    const double quad_length = 0.5;
    const double sep = 10.0;
    const double focus = 7.0;
    const double quad_strength = 1.0/(focus*quad_length);
    const double drift_length = (sep - quad_length - bend_length)/2.0;
    const double kick_strength = 0.01;
    const int ncells = 128;
    // angle = 2*pi/(2*ncells)
    const double bend_angle = mconstants::pi/ncells;
    const double co_tolerance = 1.0e-6;

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    f.set_double_attribute("k1", quad_strength);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element oh("drift", "oh");
    oh.set_double_attribute("l", 0.5*drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", -quad_strength);
    Lattice_element b("sbend", "b");
    b.set_double_attribute("l", bend_length);
    b.set_double_attribute("angle", bend_angle);

    Lattice_element k("hkicker", "k");
    k.set_double_attribute("kick", kick_strength);

    Lattice_element r("rfcavity", "r");
    r.set_double_attribute("l", 0.);
    r.set_double_attribute("harmon", 4);
    r.set_double_attribute("volt",0.0);

    Lattice_sptr lattice_sptr(new Lattice(name));
    for (int cell=0; cell<ncells; ++cell) {
        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(oh);
        lattice_sptr->append(r);
        lattice_sptr->append(oh);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
    }
    lattice_sptr->append(k);

    const int charge = pconstants::proton_charge;
    const double mass = pconstants::mp;
    const double total_energy = 2.0;
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    double beta = reference_particle.get_beta();
    lattice_sptr->set_reference_particle(reference_particle);
    double lattice_length=lattice_sptr->get_length();

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BmlPtr beamline_sptr(lattice_simulator.get_chef_lattice_sptr()->get_beamline_sptr());
    Proton probe(total_energy);
    beamline_sptr->propagate(probe);

    MArray1d closed_orbit(lattice_simulator.get_closed_orbit()); 
    MArray2d linear_map(lattice_simulator.get_linear_one_turn_map());

//     std::cout << "closed_orbit: ";
//     for (int i=0; i<6; ++i) {
//         std::cout << " " << closed_orbit[i];
//     }
//     std::cout << std::endl;
    
   // map_applied(linear_map,closed_orbit);

 //   with kick, the closed orbit in x and xp better be away from 0
    for (int i=0; i<2; ++i) {
        BOOST_CHECK(std::abs(closed_orbit[i]) > 1.0e-5);
    }

    probe.set_x(closed_orbit[0]);
    probe.set_npx(closed_orbit[1]);
    probe.set_y(closed_orbit[2]);
    probe.set_npy(closed_orbit[3]);
    probe.set_cdt(closed_orbit[4]);
    probe.set_ndp(closed_orbit[5]);
    beamline_sptr->propagate(probe);
    
   
  /*  
     std::cout << " propagate on clo before register : ";
     std::cout << probe.get_x() << " " << probe.get_npx() << " " << probe.get_y() << " " << probe.get_npy()
              << " " << probe.get_cdt() << " " << probe.get_ndp();
     std::cout << std::endl;*/

//    BOOST_CHECK(floating_point_equal(closed_orbit[0], probe.get_x(), co_tolerance));
//    BOOST_CHECK(floating_point_equal(closed_orbit[1], probe.get_npx(), co_tolerance));
//    BOOST_CHECK(floating_point_equal(closed_orbit[2], probe.get_y(), co_tolerance));
//    BOOST_CHECK(floating_point_equal(closed_orbit[3], probe.get_npy(), co_tolerance));
//    BOOST_CHECK(floating_point_equal(0., probe.get_y(), co_tolerance));
//    BOOST_CHECK(floating_point_equal(0., probe.get_npy(), co_tolerance));
   

    double off_cdt=probe.get_cdt();
    double off_clo_length=beta*off_cdt;
    double expected_rf_freq=4.*beta*pconstants::c/lattice_length;
    
   // std::cout<<" rf frequency before register="<<lattice_simulator.get_rf_frequency()<<std::endl;
   // std::cout<<" rf expected frequency before register="<<expected_rf_freq<<std::endl;
    
    BOOST_CHECK(floating_point_equal(lattice_simulator.get_rf_frequency(), expected_rf_freq, co_tolerance));

//******************************************************************************************************
//  register closed orbit    
    lattice_simulator.register_closed_orbit(); 
    double clo_length=lattice_simulator.get_closed_orbit_length();
    
    probe.set_x(closed_orbit[0]);
    probe.set_npx(closed_orbit[1]);
    probe.set_y(closed_orbit[2]);
    probe.set_npy(closed_orbit[3]);
    probe.set_cdt(closed_orbit[4]);
    probe.set_ndp(closed_orbit[5]);
    beamline_sptr->propagate(probe);
    

     MArray1d closed_orbit1(lattice_simulator.get_closed_orbit());      
     multi_array_check_equal(closed_orbit, closed_orbit1,  10.e-8);
     MArray2d linear_map1(lattice_simulator.get_linear_one_turn_map());
     multi_array_check_equal(linear_map, linear_map1,  10.e-8);

//     std::cout << " propagate on clo after register: ";
//     std::cout << probe.get_x() << " " << probe.get_npx() << " " << probe.get_y() << " " << probe.get_npy()
//               << " " << probe.get_cdt() << " " << probe.get_ndp();
//     std::cout << std::endl;
   // map_applied(linear_map1,closed_orbit1);

    BOOST_CHECK(floating_point_equal(closed_orbit[0], probe.get_x(), co_tolerance));
    BOOST_CHECK(floating_point_equal(closed_orbit[1], probe.get_npx(), co_tolerance));
    BOOST_CHECK(floating_point_equal(closed_orbit[2], probe.get_y(), co_tolerance));
    BOOST_CHECK(floating_point_equal(closed_orbit[3], probe.get_npy(), co_tolerance));
    BOOST_CHECK(floating_point_equal(0., probe.get_cdt(), co_tolerance));
    BOOST_CHECK(floating_point_equal(0., probe.get_ndp(), co_tolerance));
   
    
    expected_rf_freq=4.*beta*pconstants::c/clo_length;

//     std::ios_base::fmtflags old_flags(std::cout.flags());
//     std::cout << std::scientific << std::setprecision(8);    
//     std::cout <<"lattice_length="<<lattice_simulator.get_lattice_sptr()->get_length()<<std::endl;
//     std::cout <<"after rg closed_orbit_length="<<lattice_simulator.get_closed_orbit_length()<<std::endl;
//     std::cout <<" beta="<<beta<<std::endl;
//     std::cout <<" off_cdt="<<off_cdt<<std::endl;
//     std::cout <<" off_clo_length="<<off_clo_length<<" clo_length-lattice length="<<clo_length-lattice_length<< std::endl;
//     std::cout <<" expected_rf_freq="<<expected_rf_freq<< std::endl;
//     std::cout.flags(old_flags);
    
     BOOST_CHECK(floating_point_equal(off_clo_length, clo_length-lattice_length, co_tolerance));    
     BOOST_CHECK(floating_point_equal(lattice_simulator.get_rf_frequency(), expected_rf_freq, co_tolerance));     
     lattice_simulator.update();
     BOOST_CHECK(floating_point_equal(lattice_simulator.get_rf_frequency(), expected_rf_freq, co_tolerance));
     
     
}


BOOST_FIXTURE_TEST_CASE(calculate_element_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.calculate_element_lattice_functions();
}


BOOST_FIXTURE_TEST_CASE(calculate_slice_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*it, 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*it, 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.set_slices(slices);
    lattice_simulator.calculate_slice_lattice_functions();
}

BOOST_FIXTURE_TEST_CASE(get_element_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        Lattice_functions f(lattice_simulator.get_lattice_functions(*(*it)));
    }
}

BOOST_FIXTURE_TEST_CASE(get_slice_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*it, 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*it, 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.set_slices(slices);
    for (Lattice_element_slices::iterator it = slices.begin();
            it != slices.end(); ++it) {
        Lattice_functions f(lattice_simulator.get_lattice_functions(*(*it)));
    }
}

BOOST_FIXTURE_TEST_CASE(get_tunes, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    const double tolerance = 1.0e-3;
    const double expected_horizontal_tune = 0.70859;
    const double expected_vertical_tune = 0.00865009;
    BOOST_CHECK_CLOSE(lattice_simulator.get_horizontal_tune(),
            expected_horizontal_tune, tolerance);
    BOOST_CHECK_CLOSE(lattice_simulator.get_vertical_tune(),
            expected_vertical_tune, tolerance);
    std::pair<double, double > the_tunes = lattice_simulator.get_both_tunes();
    BOOST_CHECK_CLOSE(the_tunes.first, expected_horizontal_tune, tolerance);
    BOOST_CHECK_CLOSE(the_tunes.second, expected_vertical_tune, tolerance);
}

BOOST_FIXTURE_TEST_CASE(adjust_tunes, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Lattice_elements horizontal_correctors, vertical_correctors;
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "quadrupole") {
            if ((*it)->get_double_attribute("k1") > 0.0) {
                horizontal_correctors.push_back(*it);
            } else {
                vertical_correctors.push_back(*it);
            }
        }
    }
    const double new_horizontal_tune = 0.69;
    const double new_vertical_tune = 0.15;
    const double tolerance = 1.0e-7;
    const int verbosity = 0;
    lattice_simulator.adjust_tunes(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors, tolerance, verbosity);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_horizontal_tune() - new_horizontal_tune) < tolerance);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_vertical_tune() - new_vertical_tune) < tolerance);
}

BOOST_FIXTURE_TEST_CASE(adjust_tunes_chef, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Lattice_elements horizontal_correctors, vertical_correctors;
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "quadrupole") {
            if ((*it)->get_double_attribute("k1") > 0.0) {
                horizontal_correctors.push_back(*it);
            } else {
                vertical_correctors.push_back(*it);
            }
        }
    }
    const double new_horizontal_tune = 0.69;
    const double new_vertical_tune = 0.15;
    const double tolerance = 1.0e-7;
    const int verbosity = 0;
    const int max_steps=10;
    lattice_simulator.adjust_tunes_chef(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors, max_steps, tolerance, verbosity);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_horizontal_tune() - new_horizontal_tune) < tolerance);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_vertical_tune() - new_vertical_tune) < tolerance);
}


void print_precalc_map(MArray2d const& map)
{
    std::cout << "const double precalc_map[6][6] = {\n";
    for (int i = 0; i < 6; ++i) {
        std::cout << "    {";
        for (int j = 0; j < 6; ++j) {
            std::cout << std::setprecision(16) << map[i][j];
            if (j < 5) {
                std::cout << ", ";
            }
        }
        if (i < 5) { 
            std::cout << "},\n";
        } else {
            std::cout << "}\n";
        }
    }
    std::cout << "};\n";    
}

BOOST_FIXTURE_TEST_CASE(get_linear_one_turn_map, Foborodobo32_fixture)
{
    const int map_order = 1;
    const double tolerance = 1.0e-10;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const double precalc_map[6][6] = {
        {-2.193577261287337, 32.93854148278354, 0, 0, -5.621693373928432e-05, 2.103705558674908},
        {-0.1980015732215495, 2.517267683732674, 0, 0, -3.530199593352855e-05, 0.2250923801266718},
        {0, 0, 1.070334647703023, 1.265501306265068, 0, 0},
        {0, 0, -0.04372593897427258, 0.8825882345653975, 0, 0},
        {-0.07764401933834172, 2.126311447299118, 0, 0, 0.996935702923177, 4.907233406896752},
        {-1.786741621293265e-05, -0.0003111856575877582, 0, 0, -0.0006283185307179529, 1.000043004777226}
    };

    MArray2d calculated_map(lattice_simulator.get_linear_one_turn_map());
   // print_precalc_map(calculated_map);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            BOOST_CHECK(
                        floating_point_equal(calculated_map[i][j], precalc_map[i][j],tolerance));
        }
    }
}



BOOST_FIXTURE_TEST_CASE(get_linear_one_turn_map_after_get_tunes, Foborodobo32_fixture)
{
    const int map_order = 1;
    const double tolerance = 1.0e-10;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    // This test fails before update() is added to the get_xxxxx_tune() routine.
    const double expected_frac_tune = 0.224126196916268;
    const double expected_eigen_tune = 0.22412619691626814;

    double horizontal_tune = lattice_simulator.get_horizontal_tune();
    BOOST_CHECK_CLOSE(horizontal_tune, expected_frac_tune, tolerance);
    horizontal_tune = lattice_simulator.get_horizontal_tune(1);
    BOOST_CHECK_CLOSE(horizontal_tune, expected_eigen_tune, tolerance);


    const double precalc_map[6][6] = {
        {-2.193577261287337, 32.93854148278354, 0, 0, -5.621693373928432e-05, 2.103705558674908},
        {-0.1980015732215495, 2.517267683732674, 0, 0, -3.530199593352855e-05, 0.2250923801266718},
        {0, 0, 1.070334647703023, 1.265501306265068, 0, 0},
        {0, 0, -0.04372593897427258, 0.8825882345653975, 0, 0},
        {-0.07764401933834172, 2.126311447299118, 0, 0, 0.996935702923177, 4.907233406896752},
        {-1.786741621293265e-05, -0.0003111856575877582, 0, 0, -0.0006283185307179529, 1.000043004777226}
    };
    
    MArray2d calculated_map(lattice_simulator.get_linear_one_turn_map());
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            BOOST_CHECK(
                    floating_point_equal(calculated_map[i][j], precalc_map[i][j],tolerance));
        }
    }
}



BOOST_FIXTURE_TEST_CASE(get_chromaticities, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    const double tolerance = 1.0e-3;
    double chrH = lattice_simulator.get_horizontal_chromaticity();
    double chrV = lattice_simulator.get_vertical_chromaticity();
    const double expected_horizontal_chromaticity = -2.15053914;
    const double expected_vertical_chromaticity = -2.14414478;
    const double expected_slip_factor=0.148477934;
    const double expected_momentum_compaction=0.162233472;
    BOOST_CHECK_CLOSE(lattice_simulator.get_horizontal_chromaticity(),
            expected_horizontal_chromaticity, tolerance);
    BOOST_CHECK_CLOSE(lattice_simulator.get_vertical_chromaticity(),
            expected_vertical_chromaticity, tolerance);
    BOOST_CHECK_CLOSE(lattice_simulator.get_slip_factor(),
	    expected_slip_factor,  tolerance);  
    BOOST_CHECK_CLOSE(lattice_simulator.get_momentum_compaction(),
		 expected_momentum_compaction, tolerance);     
	    
  /*    
    std::cout<<"chromaticities (H,V):  ("<< chrH<<" ,  "<<chrV<<")"<<std::endl;
    std::cout<<" slip factor, momentum_compaction="<<lattice_simulator.get_slip_factor()
    <<"  ,  "<<lattice_simulator.get_momentum_compaction()<<std::endl;*/

}

BOOST_FIXTURE_TEST_CASE(serialize_sliced_chef_beamline, Lattice_fixture)
{
   Lattice_simulator lattice_simulator(lattice_sptr, map_order);
   Lattice_element_slices slices;
   for (Lattice_elements::const_iterator it =
           lattice_sptr->get_elements().begin(); it
           != lattice_sptr->get_elements().end(); ++it) {
       double length = (*it)->get_length();
       Lattice_element_slice_sptr first_half(
               new Lattice_element_slice(*it, 0.0, 0.5 * length));
       Lattice_element_slice_sptr second_half(
               new Lattice_element_slice(*it, 0.5 * length, length));
       slices.push_back(first_half);
       slices.push_back(second_half);
   }
   lattice_simulator.set_slices(slices);

   xml_save(lattice_simulator, "lattice_simulator2.xml");

   Lattice_simulator loaded;
   xml_load(loaded, "lattice_simulator2.xml");
}

BOOST_FIXTURE_TEST_CASE(serialize_registered_close_orbit, Lattice_fixture)
{
  
    const double bend_length = 2.0;
    const double quad_length = 0.5;
    const double sep = 10.0;
    const double focus = 7.0;
    const double quad_strength = 1.0/(focus*quad_length);
    const double drift_length = (sep - quad_length - bend_length)/2.0;
    const double kick_strength = 0.01;
    const int ncells = 128;
    // angle = 2*pi/(2*ncells)
    const double bend_angle = mconstants::pi/ncells;
    const double co_tolerance = 1.0e-6;

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    f.set_double_attribute("k1", quad_strength);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element oh("drift", "oh");
    oh.set_double_attribute("l", 0.5*drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", -quad_strength);
    Lattice_element b("sbend", "b");
    b.set_double_attribute("l", bend_length);
    b.set_double_attribute("angle", bend_angle);

    Lattice_element k("hkicker", "k");
    k.set_double_attribute("kick", kick_strength);

    Lattice_element r("rfcavity", "r");
    r.set_double_attribute("l", 0.);
    r.set_double_attribute("harmon", 4);
    r.set_double_attribute("volt",0.0);

    Lattice_sptr lattice_sptr(new Lattice(name));
    for (int cell=0; cell<ncells; ++cell) {
        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(oh);
        lattice_sptr->append(r);
        lattice_sptr->append(oh);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
    }
    lattice_sptr->append(k);
    
    const int charge = pconstants::proton_charge;
    const double mass = pconstants::mp;
    const double total_energy = 2.0;
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    double beta = reference_particle.get_beta();
    lattice_sptr->set_reference_particle(reference_particle);
    double lattice_length=lattice_sptr->get_length();
  
    
     Lattice_simulator lattice_simulator(lattice_sptr, map_order);     
     lattice_simulator.register_closed_orbit();
    

     MArray1d closed_orbit(lattice_simulator.get_closed_orbit());
     BmlPtr beamline_sptr(lattice_simulator.get_chef_lattice_sptr()->get_beamline_sptr());
     Proton probe(total_energy);


     probe.set_x(closed_orbit[0]);
     probe.set_npx(closed_orbit[1]);
     probe.set_y(closed_orbit[2]);
     probe.set_npy(closed_orbit[3]);
     probe.set_cdt(closed_orbit[4]);
     probe.set_ndp(closed_orbit[5]);
     beamline_sptr->propagate(probe);
    
 

     
     
     
     xml_save(lattice_simulator, "lattice_simulator_rclo.xml");

     Lattice_simulator loaded;
     xml_load(loaded, "lattice_simulator_rclo.xml");
     
     MArray1d closed_orbit1(loaded.get_closed_orbit());      
     multi_array_check_equal(closed_orbit, closed_orbit1,  1.e-8);     
     BOOST_CHECK(floating_point_equal(lattice_simulator.get_rf_frequency(),loaded.get_rf_frequency(),1.e-10));
     BOOST_CHECK(floating_point_equal(lattice_simulator.get_rf_bucket_length(),loaded.get_rf_bucket_length(),1.e-10));

     BmlPtr beamline1_sptr(loaded.get_chef_lattice_sptr()->get_beamline_sptr());
     Proton probe1(total_energy);


     probe1.set_x(closed_orbit[0]);
     probe1.set_npx(closed_orbit[1]);
     probe1.set_y(closed_orbit[2]);
     probe1.set_npy(closed_orbit[3]);
     probe1.set_cdt(closed_orbit[4]);
     probe1.set_ndp(closed_orbit[5]);
     beamline1_sptr->propagate(probe1);
     
     
    BOOST_CHECK(floating_point_equal(probe.get_x(), probe1.get_x(), co_tolerance));
    BOOST_CHECK(floating_point_equal(probe.get_npx(), probe1.get_npx(), co_tolerance));
    BOOST_CHECK(floating_point_equal(probe.get_y(), probe1.get_y(), co_tolerance));
    BOOST_CHECK(floating_point_equal(probe.get_npy(), probe1.get_npy(), co_tolerance));
    BOOST_CHECK(floating_point_equal(probe.get_cdt(), probe1.get_cdt(), co_tolerance));
    BOOST_CHECK(floating_point_equal(probe.get_ndp(), probe1.get_ndp(), co_tolerance)); 
   
}

BOOST_FIXTURE_TEST_CASE(is_ring, Foborodobo32_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BOOST_CHECK(lattice_simulator.is_ring());
}


BOOST_FIXTURE_TEST_CASE(adjust_chromaticities, Fosobodosobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    const double tolerance = 1.0e-3;
    double chr_h = lattice_simulator.get_horizontal_chromaticity();
    double chr_v = lattice_simulator.get_vertical_chromaticity();
//    std::cout << "begin chromaticities (H,V):  (" << chr_h << " ,  " << chr_v
//            << ")" << std::endl;

    Lattice_elements horizontal_correctors;
    Lattice_elements vertical_correctors;
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_name() == "s1") {
            // std::cout<<" h name="<<(*it)->get_name()<<" h type="<<(*it)->get_type()<<std::endl;
            horizontal_correctors.push_back(*it);
        } else if ((*it)->get_name() == "s2") {
            // std::cout<<" v name="<<(*it)->get_name()<<" v type="<<(*it)->get_type()<<std::endl;
            vertical_correctors.push_back(*it);
        }
    }
    const double newchr_h = -2.9;
    const double newchr_v = -3.1;
    const double chrom_tolerance = 5.0e-7;
    const int max_steps = 5;
    lattice_simulator.adjust_chromaticities(newchr_h, newchr_v, horizontal_correctors,
            vertical_correctors, chrom_tolerance, max_steps);

    chr_h = lattice_simulator.get_horizontal_chromaticity();
    chr_v = lattice_simulator.get_vertical_chromaticity();

//    std::cout << "final chromaticities (H,V):  (" << chr_h << " ,  " << chr_v
//            << ")" << std::endl;
    double percent_chrom_tolerance = 100 * chrom_tolerance;
    BOOST_CHECK_CLOSE(chr_h, newchr_h, percent_chrom_tolerance);
    BOOST_CHECK_CLOSE(chr_v, newchr_v, percent_chrom_tolerance);
}



BOOST_FIXTURE_TEST_CASE(xyz_normal_xyz, Foborodobo32_fixture)
{
    const int map_order = 3;
    // tolerance in transverse coordinates can be stricter than
    // longitudinal coordinates because the longitudinal coordinate
    // has the explicit sine wave which blows up when truncated.
    const double trans_tolerance = 1.0e-5;
    const double long_tolerance = 1.0e-3;

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    // fill the bunch with three points each direction
    const double test_points[] = { 1.0e-4, 7.76e-6, 1.0e-4, 1.86e-5, 1.0e-4,
            1.0e-6 };

    int pidx = 0;
    const int num_macro_particles = 3 * 3 * 3 * 3 * 3 * 3;
    MArray2d particles(boost::extents[num_macro_particles][7]);

    for (int i0 = -1; i0 != 2; ++i0) {
        for (int i1 = -1; i1 != 2; ++i1) {
            for (int i2 = -1; i2 != 2; ++i2) {
                for (int i3 = -1; i3 != 2; ++i3) {
                    for (int i4 = -1; i4 != 2; ++i4) {
                        for (int i5 = -1; i5 != 2; ++i5) {
                            particles[pidx][0] = test_points[0] * i0;
                            particles[pidx][1] = test_points[1] * i1;
                            particles[pidx][2] = test_points[2] * i2;
                            particles[pidx][3] = test_points[3] * i3;
                            particles[pidx][4] = test_points[4] * i4;
                            particles[pidx][5] = test_points[5] * i5;
                            ++pidx;
                        }
                    }
                }
            }
        }
    }

    MArray2d particles_orig(particles);

    lattice_simulator.convert_xyz_to_normal(particles);
    lattice_simulator.convert_normal_to_xyz(particles);

    for (int i = 0; i < num_macro_particles; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (j < 4) {
                BOOST_CHECK(
                        floating_point_equal(particles[i][j], particles_orig[i][j], trans_tolerance));
            } else {
                BOOST_CHECK(
                        floating_point_equal(particles[i][j], particles_orig[i][j], long_tolerance));
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(normal_xyz_normal, Foborodobo32_fixture)
{
    const int map_order = 3;
    // tolerance in transverse coordinates can be stricter than
    // longitudinal coordinates because the longitudinal coordinate
    // has the explicit sine wave which blows up when truncated.
    const double trans_tolerance = 1.0e-4;
    const double long_tolerance = 1.0e-3;

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    // fill the bunch with three points at fixed action each
    // direction but angles uniformly spread around 2*pi

    const int n_angles = 24;

    // three sets of actions
    // this is the square-root of the action
    const double action_sets[][3] = { { 6.0e-8, 0.0, 2.5e-9 }, { 0.0, 6.0e-8,
            0.0 }, { 1.9e-9, 0.0, 4.0e-8 } };

    int nsets = sizeof(action_sets) / (3 * sizeof(double));
    int pidx = 0;

    const int num_macro_particles = nsets * n_angles * n_angles * n_angles;
    MArray2d particles(boost::extents[num_macro_particles][7]);

    for (int aset = 0; aset < nsets; ++aset) {
        for (int iph0 = 0; iph0 < n_angles; ++iph0) {
            double phase0 = (2.0 * mconstants::pi / (2.0 * n_angles))
                    * (2 * iph0 + 1);
            for (int iph1 = 0; iph1 < n_angles; ++iph1) {
                double phase1 = (2.0 * mconstants::pi / (2.0 * n_angles))
                        * (2 * iph1 + 1);
                for (int iph2 = 0; iph2 < n_angles; ++iph2) {
                    double phase2 = (2.0 * mconstants::pi / (2.0 * n_angles))
                            * (2 * iph2 + 1);
                    particles[pidx][0] = action_sets[aset][0] * sin(phase0);
                    particles[pidx][1] = -action_sets[aset][0] * cos(phase0);
                    particles[pidx][2] = action_sets[aset][1] * sin(phase1);
                    particles[pidx][3] = -action_sets[aset][1] * cos(phase1);
                    particles[pidx][4] = action_sets[aset][2] * sin(phase2);
                    particles[pidx][5] = -action_sets[aset][2] * cos(phase2);
                    ++pidx;
                }
            }
        }
    }

    MArray2d particles_orig(particles);

    lattice_simulator.convert_normal_to_xyz(particles);
    lattice_simulator.convert_xyz_to_normal(particles);

    for (int i = 0; i < num_macro_particles; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (j < 4) {
                BOOST_CHECK(
                        floating_point_equal(particles[i][j], particles_orig [i][j], trans_tolerance));
            } else {
                BOOST_CHECK(
                        floating_point_equal(particles[i][j], particles_orig[i][j], long_tolerance));
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(check_linear_normal_form, Foborodobo32_fixture)
{
    const int map_order = 3;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BOOST_CHECK(lattice_simulator.check_linear_normal_form());
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    xml_save(lattice_simulator, "lattice_simulator1.xml");

//    Lattice_simulator loaded;
//    xml_load(loaded, "lattice_simulator1.xml");
}

BOOST_AUTO_TEST_CASE(lattice_functions_from_maps)
{
	const double a_x = 0.1;
	const double b_x = 32.0;
	const double mu_x = 2.0 * mconstants::pi * 0.42526;
	const double a_y = -0.8;
	const double b_y = 4.5;
	const double mu_y = 2.0 * mconstants::pi * 0.41526;
	const double b_l = 415.0;
	const double mu_l = (2.0 * mconstants::pi)/995.0;

	MArray2d map(boost::extents[6][6]);
	map[0][0] = cos(mu_x) + a_x * sin(mu_x);
	map[0][1] = b_x * sin(mu_x);
	map[1][0] = -((1 + a_x*a_x)/b_x) * sin(mu_x);
	map[1][1] = cos(mu_x) - a_x * sin(mu_x);

	map[2][2] = cos(mu_y) + a_y * sin(mu_y);
	map[2][3] = b_y * sin(mu_y);
	map[3][2] = -((1 + a_y*a_y)/b_y) * sin(mu_y);
	map[3][3] = cos(mu_y) - a_y * sin(mu_y);

	map[4][4] = cos(mu_l);
	map[4][5] = b_l * sin(mu_l);
	map[5][4] = -(1.0/b_l) * sin(mu_l);
	map[5][5] = cos(mu_l);

	Lattice_functions lf(map);
	Long_lattice_functions llf(map);

	BOOST_CHECK(floating_point_equal(a_x,lf.alpha_x, tolerance));
	BOOST_CHECK(floating_point_equal(b_x, lf.beta_x, tolerance));
	BOOST_CHECK(floating_point_equal(mu_x, lf.psi_x, tolerance));

	BOOST_CHECK(floating_point_equal(a_y,lf.alpha_y, tolerance));
	BOOST_CHECK(floating_point_equal(b_y, lf.beta_y, tolerance));
	BOOST_CHECK(floating_point_equal(mu_y, lf.psi_y, tolerance));

	BOOST_CHECK(floating_point_equal(0.0, llf.alpha, tolerance));
	BOOST_CHECK(floating_point_equal(b_l, llf.beta, tolerance));
	BOOST_CHECK(floating_point_equal(mu_l, llf.psi, tolerance));

}


