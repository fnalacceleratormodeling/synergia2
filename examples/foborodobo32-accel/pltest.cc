
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/populate.h"
#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/foundation/pcg_distribution.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/simulation/checkpoint.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/populate_stationary.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/utils/lsexpr.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/utils.h"

using LV = LoggerV;

Lattice
get_lattice()
{
    namespace LS = Lattice_simulator;

    MadX_reader reader;
    Lattice lattice(reader.get_lattice("model", "cfoborodobo32.madx"));

    lattice.set_all_string_attribute("extractor_type", "libff");

    Reference_particle ref(1, pconstants::mp, pconstants::mp + 0.8);

    lattice.set_reference_particle(ref);

    // tune the lattice
    LS::tune_circular_lattice(lattice);

    return lattice;
}

void
run()
{

    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_STEP);

    Lattice lattice(get_lattice());

    // get the reference particle
    auto const& ref = lattice.get_reference_particle();

    screen(LV::INFO) << "reference energy: " << ref.get_total_energy();
    screen(LV::INFO) << "reference momentum: " << ref.get_momentum()
                     << " GeV\n";
    screen(LV::INFO) << "gamma: " << ref.get_gamma() << "\n";
    screen(LV::INFO) << "beta: " << ref.get_beta() << "\n";

    // stepper
    Split_operator_stepper_elements stepper(Dummy_CO_options(), 1);

    // Propagator
    Propagator propagator(lattice, stepper);

    screen(LV::INFO)
        << "lattice.get_elements() first element orig a1 parameter: "
        << lattice.get_elements().front().get_double_attribute("a1")
        << std::endl;
    lattice.get_elements().front().set_double_attribute("a1", 9999.0);
    screen(LV::INFO)
        << "lattice.get_elements() first element orig a1 parameter: "
        << lattice.get_elements().front().get_double_attribute("a1")
        << std::endl;
    screen(LV::INFO) << "demonstrated that I can control lattice element "
                        "attributes from lattice.get_elements()\n"
                     << std::endl;

    screen(LV::INFO)
        << "Initial a1 param propagator.get_lattice_elements().front(): "
        << propagator.get_lattice_elements().front().get_double_attribute("a1")
        << std::endl;
    screen(LV::INFO) << "propagator.get_lattice_elements() returns a const "
                        "reference to lattice elements that are different than "
                        "the original lattice used in propagator.  You can't "
                        "modify any elements through this reference.\n"
                     << std::endl;

    // Doesn't work because of const
    // propagator.get_lattice_elements().front().set_double_attribute("a1",
    // -888.0); screen(LV::INFO) << "After mod a1 param
    // propagator.get_lattice_elements().front()" <<
    //     propagator.get_lattice_elements().front().get_double_attribute("a1")
    //     << std::endl;

    Propagator::Lattice_element_slices les =
        propagator.get_lattice_element_slices();
    auto& first_slice = *les.begin();

    // By the way, Lattice_element_slice::get_lattice_element is a const& so
    // can't modify element?
    screen(LV::INFO) << "First slice:\n"
                     << first_slice.get_lattice_element().get_name() << ": a1: "
                     << first_slice.get_lattice_element().get_double_attribute(
                            "a1")
                     << std::endl;

    // Doesn't work because of const
    // first_slice.get_lattice_element().set_double_attribute("a1", -222.0);

    Lattice& l(propagator.get_lattice());
    screen(LV::INFO)
        << "propagator.get_lattice().get_elements().front() initial a1 param "
           "l.get_elements().front(): "
        << propagator.get_lattice().get_elements().front().get_double_attribute(
               "a1")
        << std::endl;
    propagator.get_lattice().get_elements().front().set_double_attribute(
        "a1", -555.0);
    screen(LV::INFO)
        << "propagator.get_lattice().get_elements().front() after mod a1 param "
           "l.get_elements().front(): "
        << propagator.get_lattice().get_elements().front().get_double_attribute(
               "a1")
        << std::endl;

    screen(LV::INFO) << "First slice after mod:\n"
                     << first_slice.get_lattice_element().get_name() << ": a1: "
                     << first_slice.get_lattice_element().get_double_attribute(
                            "a1")
                     << std::endl;

    return;
}

int
main(int argc, char** argv)
{
    synergia::initialize(argc, argv);

    run();

    synergia::finalize();
    return 0;
}
