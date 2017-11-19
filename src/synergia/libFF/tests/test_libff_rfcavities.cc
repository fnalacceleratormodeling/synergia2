#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/floating_point.h"
#include "synergia/utils/serialization.h"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/lattice/madx_adaptor_map.h"
#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include "synergia/bunch/bunch.h"

#include "synergia/lattice/madx_reader.h"

#include "synergia/utils/boost_test_mpi_fixture.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture); // needed to initialize MPI

#if 0
void slice_times(Lattice_simulator &lattice_simulator, fstream &fs);
#endif

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(tune_circular_lattice_cavities)
{
    MadX_reader reader;
    Lattice_sptr lattice_sptr(reader.get_lattice_sptr("model", "lattices/foborodobo128.madx"));
    lattice_sptr->set_all_string_attribute("extractor_type", "libff");
    Reference_particle reference_particle(lattice_sptr->get_reference_particle());

#if 0
    std::fstream ofs("barlattice.txt", std::fstream::out);
    ofs << lattice_sptr->as_string() << std::endl;
    ofs.close();
#endif

    const double beta = reference_particle.get_beta();
    const double momentum = reference_particle.get_momentum();
    const double lattice_length = lattice_sptr->get_length();

    // std::cout << "tune_circular_lattice, length: " << lattice_length << ", momentum: " << momentum << ", beta: " << beta << std::endl;

    // harmonic number of test lattice is 1
    const double freq = 1.0 * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    // std::cout << "before tune circular" << std::endl;
    MArray1d costate = stepper_sptr->get_lattice_simulator().tune_circular_lattice();

#if 0
    fstream sts("slice_times1.txt", std::fstream::out);
    slice_times(stepper_sptr->get_lattice_simulator(), sts);
    sts.close();
    std::cout << "closed orbit state: " << std::setprecision(14) << costate[0] << " " << costate[1] << " " << costate[2] << " " << costate[3] << " " << costate[4] << " " << costate[5] << std::endl;
    std::cout << "after tune circular" << std::endl;
#endif
    // stepper_sptr->get_lattice_simulator().tune_circular_lattice();
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_linear_lattice_cavities)
{
    const double length = 10.0;
    const double rfclength = 0.5;
    const double rfcvolt = 0.01;
    const double rfcharmon = 10.0;
    const double pc = 0.5;
    
    // create lattice
    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element rfc("rfcavity", "rfc");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    lattice_sptr->append(rfc);
    
    Lattice_element d("drift", "food");
    d.set_double_attribute("l", length);
    lattice_sptr->append(d);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    const double beta = reference_particle.get_beta();
    const double lattice_length = lattice_sptr->get_length();

    const double freq = rfcharmon * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    stepper_sptr->get_lattice_simulator().tune_linear_lattice();
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_diagonal_lattice_cavities)
{
    const double length = 10.0;
    const double rfclength = 0.5;
    const double rfcvolt = 0.01;
    const double rfcharmon = 10.0;
    const double pc = 0.5;
    
    const double toffs = 0.01;
    
    // create lattice
    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element rfc("rfcavity", "rfc");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    lattice_sptr->append(rfc);
    
    Lattice_element d("drift", "food");
    d.set_double_attribute("l", length);
    lattice_sptr->append(d);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    
    // set momentum for diagonal path
    const double lattice_length = lattice_sptr->get_length();
    double total_length = std::sqrt(lattice_length*lattice_length + toffs*toffs);
    MArray1d state = reference_particle.get_state();
    state[Bunch::yp] = toffs/total_length;
    reference_particle.set_state(state);
    
    lattice_sptr->set_reference_particle(reference_particle);

    const double beta = reference_particle.get_beta();

    const double freq = rfcharmon * beta * pconstants::c/total_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    stepper_sptr->get_lattice_simulator().tune_linear_lattice();
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_bend_lattice_cavities)
{
    const double length = 10.0;
    const double angle = 0.01;
    const double rfclength = 0.5;
    const double rfcvolt = 0.01;
    const double rfcharmon = 10.0;
    const double pc = 0.5;
    
    // create lattice
    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element rfc("rfcavity", "rfc");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    lattice_sptr->append(rfc);
    
    Lattice_element b("sbend", "foos");
    b.set_double_attribute("l", length);
    b.set_double_attribute("angle", angle);
    lattice_sptr->append(b);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    const double beta = reference_particle.get_beta();
    const double lattice_length = lattice_sptr->get_length();

    const double freq = rfcharmon * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    stepper_sptr->get_lattice_simulator().tune_linear_lattice();
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_quad_lattice_cavities)
{
    const double length = 10.0;
    const double rfclength = 0.5;
    const double rfcvolt = 0.01;
    const double rfcharmon = 10.0;
    const double qlength = 0.5;
    const double qstrength = 0.2;
    const double pc = 0.5;
    
    // create lattice
    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element rfc("rfcavity", "rfc");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    lattice_sptr->append(rfc);

    Lattice_element q("quadrupole", "q");
    q.set_double_attribute("l", qlength);
    q.set_double_attribute("k1", qstrength);
    lattice_sptr->append(q);
    
    Lattice_element d("drift", "food");
    d.set_double_attribute("l", length);
    lattice_sptr->append(d);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    const double beta = reference_particle.get_beta();
    const double lattice_length = lattice_sptr->get_length();

    const double freq = rfcharmon * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    stepper_sptr->get_lattice_simulator().tune_linear_lattice();
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_circular_drift_cavities)
{
    const double length = 2560.0;
    const double rfclength = 0.5;
    const double rfcvolt = 0.01;
    const double rfcharmon = 10.0;
    const double pc = 0.5;
    
    MadX_adaptor_map_sptr madx_adaptor_map_sptr(new MadX_adaptor_map);

    // create lattice
    Lattice_sptr lattice_sptr(new Lattice("foo", madx_adaptor_map_sptr));

    Lattice_element rfc("rfcavity", "rfc");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    lattice_sptr->append(rfc);
    
    Lattice_element d("drift", "food");
    d.set_double_attribute("l", length);
    lattice_sptr->append(d);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    const double beta = reference_particle.get_beta();
    const double lattice_length = lattice_sptr->get_length();

    const double freq = rfcharmon * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    std::cout << "before tune circular" << std::endl;
    MArray1d costate = stepper_sptr->get_lattice_simulator().tune_circular_lattice();
    std::cout << "closed orbit state: " << std::setprecision(14) << costate[0] << " " << costate[1] << " " << costate[2] << " " << costate[3] << " " << costate[4] << " " << costate[5] << std::endl;
    std::cout << "after tune circular" << std::endl;
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_circular_monster_cavities)
{
    const double length = 10.0;
    const double rfclength = 0.5;
    const double rfcvolt = 0.01;
    const double rfcharmon = 10.0;
    const double qlength = 0.5;
    const double qstrength = 0.2;
    const double pc = 0.5;
    
    // create lattice
    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element rfc("rfcavity", "rfc");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    lattice_sptr->append(rfc);
    
    Lattice_element hk("hkicker", "hk");
    hk.set_double_attribute("l", 0.0);
    hk.set_double_attribute("kick", 0.0);
    lattice_sptr->append(hk);
    
    Lattice_element d("drift", "food");
    d.set_double_attribute("l", length);
    lattice_sptr->append(d);

    Lattice_element vk("vkicker", "vk");
    vk.set_double_attribute("l", 0.0);
    vk.set_double_attribute("kick", 0.0);
    lattice_sptr->append(vk);

    Lattice_element q("quadrupole", "q");
    q.set_double_attribute("l", qlength);
    q.set_double_attribute("k1", qstrength);
    lattice_sptr->append(q);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    const double beta = reference_particle.get_beta();
    const double lattice_length = lattice_sptr->get_length();

    const double freq = rfcharmon * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    std::cout << "before tune circular" << std::endl;
    MArray1d costate = stepper_sptr->get_lattice_simulator().tune_linear_lattice();
    std::cout << "closed orbit state: " << std::setprecision(14) << costate[0] << " " << costate[1] << " " << costate[2] << " " << costate[3] << " " << costate[4] << " " << costate[5] << std::endl;
    std::cout << "after tune circular" << std::endl;
    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(tune_circular_foborodobo128_cavities)
{
    // creating foborodobo128 lattice by hand
    const int ncell = 128;
    const double sepn = 10.0;
    const double focus = 7.0;
    const double quadlength = 2.0;
    const double strength = 1.0/(focus*quadlength);
    const double pct = 0.5;
    const double bendlength = pct*(sepn-quadlength);
    const double driftlength = (sepn-quadlength-bendlength)/2;
    const double bendangle = 2.0*mconstants::pi/(2*ncell);
    
    const double rfclength = 0.0;
    const double rfcvolt = 1.0e-4;
    const double rfcharmon = 1.0;
    
    const double pc = 0.5;
    
    //MadX_adaptor_map_sptr madx_adaptor_map_sptr;
    // create lattice
    //Lattice_sptr lattice_sptr(new Lattice("foo", madx_adaptor_map_sptr));
    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element rfc("rfcavity", "r");
    rfc.set_double_attribute("l", rfclength);
    rfc.set_double_attribute("volt", rfcvolt);
    rfc.set_double_attribute("harmon", rfcharmon);
    rfc.set_double_attribute("lag", 0.0);

    Lattice_element b("sbend", "b");
    b.set_double_attribute("l", bendlength);
    b.set_double_attribute("angle", bendangle);
    
    Lattice_element d("drift", "o");
    d.set_double_attribute("l", driftlength);

    Lattice_element shortd("drift", "os");
    shortd.set_double_attribute("l", driftlength/2.0);

    Lattice_element fq("quadrupole", "f");
    fq.set_double_attribute("l", quadlength);
    fq.set_double_attribute("k1", strength);

    Lattice_element dq("quadrupole", "d");
    dq.set_double_attribute("l", quadlength);
    dq.set_double_attribute("k1", -strength);
    
    Lattice_element vk("vkicker", "vk");
    vk.set_double_attribute("l", 0.0);
    vk.set_double_attribute("kicker", 0.0);
    
    Lattice_element hk("hkicker", "hk");
    hk.set_double_attribute("l", 0.0);
    hk.set_double_attribute("kicker", 0.0);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    // build up the full lattice as fobodobo the first cell having an
    // rf cavity in the middle of the second drift
    for (int cell=0; cell<ncell; ++cell) {
        stringstream ss;
        ss << "cell_";
        ss << std::setw(3) << std::setfill('0') << cell;
        Lattice_element le_cell("marker", ss.str());
        lattice_sptr->append(le_cell);
        lattice_sptr->append( (cell%2 == 0) ? hk : vk);
        lattice_sptr->append(fq);
        lattice_sptr->append(d);
        lattice_sptr->append(b);
        if (cell == 0) {
            lattice_sptr->append(shortd);
            lattice_sptr->append(rfc);
            lattice_sptr->append(shortd);
        } else {
            lattice_sptr->append(d);
        }
        lattice_sptr->append(dq);
        lattice_sptr->append(d);
        lattice_sptr->append(b);
        lattice_sptr->append(d);
    }
    
    lattice_sptr->set_all_string_attribute("extractor_type", "libff");

#if 0
    std::fstream ofs("foolattice.txt", std::fstream::out);
    ofs << lattice_sptr->as_string() << std::endl;
    ofs.close();
#endif

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    const double momentum = reference_particle.get_momentum();
    const double beta = reference_particle.get_beta();
    const double lattice_length = lattice_sptr->get_length();

    // std::cout << "tune_circular_foborodobo128, length: " << lattice_length << ", momentum: " << momentum << ", beta: " << beta << std::endl;
    const double freq = rfcharmon * beta * pconstants::c/lattice_length;
    
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));
    // tune the cavities
    // std::cout << "before tune circular" << std::endl;
    MArray1d costate = stepper_sptr->get_lattice_simulator().tune_circular_lattice();
    // std::cout << "closed orbit state: " << std::setprecision(14) << costate[0] << " " << costate[1] << " " << costate[2] << " " << costate[3] << " " << costate[4] << " " << costate[5] << std::endl;
    // std::cout << "after tune circular" << std::endl;
    
#if 0
    fstream sts("slice_times2.txt", std::fstream::out);
    slice_times(stepper_sptr->get_lattice_simulator(), sts);
    sts.close();
#endif

    // check that the frequencies of the cavities are correct
    Lattice_elements&  elements = lattice_sptr->get_elements();
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        if ((*lit)->get_type() == "rfcavity") {
            BOOST_CHECK_CLOSE((*lit)->get_double_attribute("freq"), freq*1.0e-6, tolerance);
        }
    }
}

#if 0
BOOST_AUTO_TEST_CASE(test_wtf)
{
    MadX_reader reader;
    Lattice_sptr lattice_sptr(reader.get_lattice_sptr("machine", "testlattice.madx"));
    Reference_particle reference_particle(lattice_sptr->get_reference_particle());
    
    const double beta = reference_particle.get_beta();
    std::cout << "test_wtf, beta: " << beta << std::endl;
}
#endif

#if 0
void slice_times(Lattice_simulator &lattice_simulator, fstream &fs)
{

    Lattice_element_slices::const_iterator sit = lattice_simulator.get_slices().begin();
    for (; sit != lattice_simulator.get_slices().end(); ++sit) {
        fs << (*sit)->get_lattice_element().get_name() << ", " ;
        fs << (*sit)->get_lattice_element().get_type() << ", " ;
        fs << (*sit)->get_reference_ct() << std::endl;
    }
}
#endif
