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


#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "synergia/libFF/ff_element.h"

namespace
{
    struct Closed_orbit_params
    {
        const double dpp;
        Lattice lattice;

        Closed_orbit_params(double dpp, Lattice lattice) 
            : dpp(dpp), lattice(lattice)
        {
            // turn off any RF cavities because they 
            // screw up the closed orbit calcation
            for (auto & ele : lattice.get_elements())
                if (ele.get_type() == element_type::rfcavity)
                    ele.set_double_attribute("volt", 0.0);
        }
    };

    int
    propagate_co_try(const gsl_vector *co_try, void *params, gsl_vector *co_results)
    {
        Closed_orbit_params *copp = static_cast<Closed_orbit_params *>(params);

        auto comm = Commxx();
        Bunch bunch(copp->lattice.get_reference_particle(), 
                comm.size(), 1.0e10, comm);

        auto lp = bunch.get_host_particles();

        // set phase space coordinates in bunch
        double orbit_start[4];
        orbit_start[0] = lp(0, 0) = gsl_vector_get(co_try, 0);
        orbit_start[1] = lp(0, 1) = gsl_vector_get(co_try, 1);
        orbit_start[2] = lp(0, 2) = gsl_vector_get(co_try, 2);
        orbit_start[3] = lp(0, 3) = gsl_vector_get(co_try, 3);

        lp(0, 4) = 0.0;
        lp(0, 5) = copp->dpp;

        // checkin
        bunch.checkin_particles();
        
        // propagate
        for(auto const& ele : copp->lattice.get_elements())
            FF_element::apply(ele, bunch);

        // checkout
        bunch.checkin_particles();

        gsl_vector_set(co_results, 0, lp(0, 0) - orbit_start[0]);
        gsl_vector_set(co_results, 1, lp(0, 1) - orbit_start[1]);
        gsl_vector_set(co_results, 2, lp(0, 2) - orbit_start[2]);
        gsl_vector_set(co_results, 3, lp(0, 3) - orbit_start[3]);

        return GSL_SUCCESS;
    }
}

std::array<double, 6> 
Lattice_simulator::calculate_closed_orbit(Lattice const& lattice,
        double dpp, double tolerance)
{
    // create params object, make a copy of the lattice
    Closed_orbit_params cop(dpp, lattice);

    const size_t ndim = 4; // solve closed orbit in x, xp, y, yp

    //const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid;
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    //const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_dnewton;
    gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(T, ndim);

    // co_try are the coordinates of the closed orbit
    gsl_vector *co_try = gsl_vector_alloc(ndim);

    // initialize the closed orbit
    for (int i=0; i<ndim; ++i)
        gsl_vector_set(co_try, i, 0.0);

    gsl_multiroot_function F;
    F.f = &propagate_co_try;
    F.n = ndim;
    F.params = &cop;

    gsl_multiroot_fsolver_set(solver, &F, co_try);

    int niter=0;
    const int maxiter = 100;
    gsl_vector *fvalues;

    do {
        int rc;
        rc = gsl_multiroot_fsolver_iterate(solver);

        switch(rc) 
        {
        case GSL_ENOPROG:
            throw std::runtime_error(
                    "Closed orbit solver unable to converge.  Is the tolerance too tight?");
            break;

        case GSL_EBADFUNC:
            throw std::runtime_error(
                    "Closed orbit solver failed to evaluate solution");
            break;

        default:
            break;
        }

        fvalues = gsl_multiroot_fsolver_f(solver);

    } while ((gsl_multiroot_test_residual(solver->f, tolerance) == GSL_CONTINUE) 
            && (++niter < maxiter));

    if (niter == maxiter) 
    {
        std::stringstream sstr;
        sstr << "Could not locate closed orbit after " << maxiter << " iterations";
        throw std::runtime_error( sstr.str() );
    }

    std::array<double, 6> costate;
    gsl_vector *froots = gsl_multiroot_fsolver_root(solver);
    costate[0] = gsl_vector_get(froots, 0);
    costate[1] = gsl_vector_get(froots, 1);
    costate[2] = gsl_vector_get(froots, 2);
    costate[3] = gsl_vector_get(froots, 3);
    costate[4] = 0.0;
    costate[5] = cop.dpp;

    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(co_try);
    return costate;
}



