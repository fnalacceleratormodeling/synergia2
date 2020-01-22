#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/independent_stepper_elements.h"



std::array<double, 6>
Lattice_simulator::tune_linear_lattice(Lattice & lattice)
{
    return tune_rfcavities(lattice);
}

std::array<double, 6>
Lattice_simulator::tune_circular_lattice(Lattice & lattice, double tolerance)
{
    // calculate closed orbit
    //auto state = calculate_closed_orbit(lattice, 0.0, tolerance);
    //lattice.get_reference_particle().set_state(state);
    
    return tune_rfcavities(lattice);
}

std::array<double, 6>
Lattice_simulator::tune_rfcavities(Lattice & lattice)
{
    // make a copy of the original lattice
    Lattice temp_lattice(lattice);
    auto & ref = temp_lattice.get_reference_particle();

    // set rfcavity volt to 0 on the copied lattice
    for (auto & ele : temp_lattice.get_elements())
    {
        if (ele.get_type() == element_type::rfcavity) 
            ele.set_double_attribute("volt", 0.0);
    }

    // setup the propagator
    Propagator propagator(temp_lattice, Independent_stepper_elements(1));

    // bunch simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(ref, 1, 1e09, Commxx());

    // propagate actions
    double accum_cdt = 0.0;
    sim.reg_prop_action_step_end([&accum_cdt](Bunch_simulator& sim, Lattice&, int, int, void*) { 
        accum_cdt += sim.get_bunch().get_design_reference_particle().get_state()[4]; 
    }, nullptr );

    // propagate to get the accumulated cdt
    Logger simlog(0, LoggerV::ERROR);
    propagator.propagate(sim, simlog, 1);

    // return the state
    auto state = sim.get_bunch().get_reference_particle().get_state();
    state[Bunch::cdt] = accum_cdt;

    // go back and set the frequency of the cavities based on the accumulated cdt
    double f = pconstants::c/accum_cdt;

    for (auto & ele : lattice.get_elements())
    {
        if (ele.get_type() == element_type::rfcavity) 
        {
            // set the frequency of the cavity if it doesn't already have one set and
            // there is a reasonable harmonic number
            if ( ele.get_double_attribute("freq", -1.0) <= 0.0 
                    && ele.get_double_attribute("harmon", -1.0) > 0.0 ) 
            {
                double harmon = ele.get_double_attribute("harmon");
                // MAD-X definition of frequency is MHz
                ele.set_double_attribute("freq", harmon*f*1.0e-6);
            }
        }
    }

    return state;
}



