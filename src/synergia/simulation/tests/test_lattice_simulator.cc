#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/lattice_simulator.h"
#include "lattice_fixture.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;
const int map_order = 2;

#if 0
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

#endif

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
    const double tolerance = 1.0e-6;
    lattice_simulator.adjust_tunes(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors, tolerance);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_horizontal_tune() - new_horizontal_tune) < tolerance);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_vertical_tune() - new_vertical_tune) < tolerance);
}

BOOST_FIXTURE_TEST_CASE(adjust_tunes_jfa, Fobodobo_sbend_fixture)
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
    const double tolerance = 1.0e-6;
    lattice_simulator.adjust_tunes_jfa(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors, tolerance);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_horizontal_tune() - new_horizontal_tune) < tolerance);
    BOOST_CHECK(
            std::abs(lattice_simulator.get_vertical_tune() - new_vertical_tune) < tolerance);
}

#if 0
BOOST_FIXTURE_TEST_CASE(get_linear_one_turn_map, Foborodobo32_fixture)
{
    const int map_order = 1;
    const double tolerance = 1.0e-10;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const double precalc_map[6][6] = { { -2.19357726128732, 32.9385414827834, 0,
            0, -5.62169337392918e-05, 2.1037055586748 }, { -0.198001573221548,
            2.51726768373267, 0, 0, -3.53019959335299e-05, 0.225092380126584 },
            { 0, 0, 1.07033464770303, 1.26550130626506, 0, 0 }, { 0, 0,
                    -0.043725938974272, 0.882588234565397, 0, 0 }, {
                    -0.077644019330161, 2.12631144692458, 0, 0,
                    0.996935702805962, 4.9072335958152 }, {
                    -1.78674162102745e-05, -0.000311185657541453, 0, 0,
                    -0.000628318530717954, 1.00004300477563 } };

    MArray2d gotten_map(lattice_simulator.get_linear_one_turn_map());
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            BOOST_CHECK(
                    floating_point_equal(gotten_map[i][j], precalc_map[i][j],tolerance));
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
    const double expected_eigen_tune = 0.224113175247965;

    double horizontal_tune = lattice_simulator.get_horizontal_tune();
    BOOST_CHECK_CLOSE(horizontal_tune, expected_frac_tune, tolerance);
    horizontal_tune = lattice_simulator.get_horizontal_tune(1);
    BOOST_CHECK_CLOSE(horizontal_tune, expected_eigen_tune, tolerance);


    const double precalc_map[6][6] = { { -2.19357726128732, 32.9385414827834, 0,
            0, -5.62169337392918e-05, 2.1037055586748 }, { -0.198001573221548,
            2.51726768373267, 0, 0, -3.53019959335299e-05, 0.225092380126584 },
            { 0, 0, 1.07033464770303, 1.26550130626506, 0, 0 }, { 0, 0,
                    -0.043725938974272, 0.882588234565397, 0, 0 }, {
                    -0.077644019330161, 2.12631144692458, 0, 0,
                    0.996935702805962, 4.9072335958152 }, {
                    -1.78674162102745e-05, -0.000311185657541453, 0, 0,
                    -0.000628318530717954, 1.00004300477563 } };

    MArray2d gotten_map(lattice_simulator.get_linear_one_turn_map());
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            BOOST_CHECK(
                    floating_point_equal(gotten_map[i][j], precalc_map[i][j],tolerance));
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
    const double expected_horizontal_chromaticity = -2.1483864;
    const double expected_vertical_chromaticity = -2.142446;
    BOOST_CHECK_CLOSE(lattice_simulator.get_horizontal_chromaticity(),
            expected_horizontal_chromaticity, tolerance);
    BOOST_CHECK_CLOSE(lattice_simulator.get_vertical_chromaticity(),
            expected_vertical_chromaticity, tolerance);
//     std::cout<<"chromaticities (H,V):  ("<< chrH<<" ,  "<<chrV<<")"<<std::endl;

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
    lattice_simulator.adjust_chromaticities(newchr_h, newchr_v, horizontal_correctors,
            vertical_correctors, 1.0e-6, 5);

    chr_h = lattice_simulator.get_horizontal_chromaticity();
    chr_v = lattice_simulator.get_vertical_chromaticity();

//    std::cout << "final chromaticities (H,V):  (" << chr_h << " ,  " << chr_v
//            << ")" << std::endl;
    const double chrom_tolerance = 5.0e-7;
    BOOST_CHECK_CLOSE(chr_h, newchr_h, chrom_tolerance);
    BOOST_CHECK_CLOSE(chr_v, newchr_v, chrom_tolerance);
}

BOOST_FIXTURE_TEST_CASE(is_ring, Foborodobo32_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BOOST_CHECK(lattice_simulator.is_ring());
}

BOOST_FIXTURE_TEST_CASE(human_normal_human, Foborodobo32_fixture)
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

    lattice_simulator.convert_human_to_normal(particles);
    lattice_simulator.convert_normal_to_human(particles);

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

BOOST_FIXTURE_TEST_CASE(normal_human_normal, Foborodobo32_fixture)
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

    lattice_simulator.convert_normal_to_human(particles);
    lattice_simulator.convert_human_to_normal(particles);

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

//BOOST_FIXTURE_TEST_CASE(serialize_sliced_chef_beamline, Lattice_fixture)
//{
//    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
//    Lattice_element_slices slices;
//    for (Lattice_elements::const_iterator it =
//            lattice_sptr->get_elements().begin(); it
//            != lattice_sptr->get_elements().end(); ++it) {
//        double length = (*it)->get_length();
//        Lattice_element_slice_sptr first_half(
//                new Lattice_element_slice(*it, 0.0, 0.5 * length));
//        Lattice_element_slice_sptr second_half(
//                new Lattice_element_slice(*it, 0.5 * length, length));
//        slices.push_back(first_half);
//        slices.push_back(second_half);
//    }
//    lattice_simulator.set_slices(slices);
//
//    xml_save(lattice_simulator, "lattice_simulator2.xml");
//
//    Lattice_simulator loaded;
//    xml_load(loaded, "lattice_simulator2.xml");
//}
#endif
