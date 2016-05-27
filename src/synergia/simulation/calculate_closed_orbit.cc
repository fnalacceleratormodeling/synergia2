#include <list>
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/commxx.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/propagator.h"

#include "calculate_closed_orbit.h"

#include <stdexcept>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


#define DEBUG 1

struct Closed_orbit_params
{
    const double dpp;
    Lattice_sptr working_lattice_sptr;
    Independent_stepper_sptr stepper_sptr;
    Propagator propagator;
    Commxx_sptr commxx_sptr;
    Bunch_sptr bunch_sptr;
#if 1
    MArray2d_ref lp;
#endif
    Bunch_simulator bunch_simulator;
    Closed_orbit_params(const double dpp, Lattice_sptr lattice_sptr) :
        dpp(dpp)
        , working_lattice_sptr(new Lattice(*lattice_sptr))
  #if 1
        , stepper_sptr(new Independent_stepper(working_lattice_sptr, 1, 1))
        , propagator(stepper_sptr)
        , commxx_sptr(new Commxx())
        , bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(),
                1, 1.0e10, commxx_sptr))
        , lp(bunch_sptr->get_local_particles())
        , bunch_simulator(bunch_sptr)
  #endif
    {
#if DEBUG
        std::cout << "egs: Closed_orbit_params constructor, dpp" << dpp << std::endl;
        std::cout << "egs: lattice length: " << working_lattice_sptr->get_length() << std::endl;
#endif

        lp[0][Bunch::x] = 0.0;
        lp[0][Bunch::xp] = 0.0;
        lp[0][Bunch::y] = 0.0;
        lp[0][Bunch::yp] = 0.0;
        lp[0][Bunch::cdt] = 0.0;
        lp[0][Bunch::dpop] = dpp;

        Lattice_elements elements(working_lattice_sptr->get_elements());

        // turn off any RF cavities because they screw up the closed orbit calcation
        for (Lattice_elements::const_iterator lep=elements.begin(); lep!=elements.end(); ++lep) {
            if ((*lep)->get_name() == "rfcavity") {
                (*lep)->set_double_attribute("volt", 0.0);
            }
        }
    }
};

int
propagate_co_try(const gsl_vector *co_try, void *params, gsl_vector *co_results)
{
    // set phase space coordinates in bunch
    Closed_orbit_params *copp = static_cast<Closed_orbit_params *>(params);
    double orbit_start[4];
    orbit_start[0] = copp->lp[0][0] = gsl_vector_get(co_try, 0);
    orbit_start[1] = copp->lp[0][1] = gsl_vector_get(co_try, 1);
    orbit_start[2] = copp->lp[0][2] = gsl_vector_get(co_try, 2);
    orbit_start[3] = copp->lp[0][3] = gsl_vector_get(co_try, 3);
#if DEBUG
    std::cout << "egs: propagate_co_try: initial: " << copp->lp[0][0]<<", "<<copp->lp[0][1]<<", "<<copp->lp[0][2]<<", "<<copp->lp[0][3]<<std::endl;
#endif
    copp->propagator.propagate(copp->bunch_simulator, 1, 1, DEBUG);
#if DEBUG
    std::cout << "egs: propagate_co_try: final: " << copp->lp[0][0]<<", "<<copp->lp[0][1]<<", "<<copp->lp[0][2]<<", "<<copp->lp[0][3]<<std::endl;
#endif
    gsl_vector_set(co_results, 0, std::abs(copp->lp[0][0] - orbit_start[0]));
    gsl_vector_set(co_results, 1, std::abs(copp->lp[0][1] - orbit_start[1]));
    gsl_vector_set(co_results, 2, std::abs(copp->lp[0][2] - orbit_start[2]));
    gsl_vector_set(co_results, 3, std::abs(copp->lp[0][3] - orbit_start[3]));
    return GSL_SUCCESS;
}

MArray1d
calculate_closed_orbit(const Lattice_sptr lattice_sptr, const double dpp, const double tolerance)
{
#if DEBUG
    std::cout << "egs: calculate_closed_orbit, tolerance=" << tolerance << std::endl;
    std::cout << "egs: using lattice: " << lattice_sptr->get_name() << " length: " << lattice_sptr->get_length() << ", number elements: " << lattice_sptr->get_elements().size() << std::endl;
    std::cout << "egs: beam energy: " << lattice_sptr->get_reference_particle().get_total_energy() << std::endl;
#endif

    Closed_orbit_params cop(dpp, lattice_sptr);
#if DEBUG
    std::cout << "egs: after instantiating Closed_orbit_params cop" << std::endl;
#endif
    const size_t ndim = 4; // solve closed orbit in x, xp, y, yp
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(T, ndim);
#if DEBUG
    std::cout << "egs: after gsl_multiroot_fsolver_alloc" << std::endl;
#endif
#if DEBUG
    gsl_vector *co_try = gsl_vector_alloc(ndim);
#endif
    gsl_multiroot_function F;
    F.f = &propagate_co_try;
    F.n = ndim;
    F.params = &cop;

    gsl_multiroot_fsolver_set(solver, &F, co_try);
#if DEBUG
    std::cout << "after gsl_multiroot_fsolver_set" << std::endl;
#endif
    int niter=0;
    const int maxiter = 100;
    gsl_vector *fvalues;
    do {
        int rc;
#if DEBUG
        std::cout << "egs: starting iter " << niter << std::endl;
#endif
        rc = gsl_multiroot_fsolver_iterate(solver);
        switch(rc) {
        case GSL_ENOPROG:
#if DEBUG
            std::cout << "egs: ENOPROG" << std::endl;
#endif
            throw std::runtime_error("ENOPROG");
            break;
        case GSL_EBADFUNC:
#if DEBUG
            std::cout << "egs: EBADFUNC" << std::endl;
#endif
            throw std::runtime_error("EBADFUNC");
            break;
        default:
#if DEBUG
            std::cout << "egs: success" << std::endl;
#endif
            break;
        }
        fvalues = gsl_multiroot_fsolver_f(solver);
#if DEBUG
        std::cout << "egs: values at current iteration: " << gsl_vector_get(fvalues, 0) << ", "
            << gsl_vector_get(fvalues, 1) << ", " << gsl_vector_get(fvalues, 2) <<
            ", " << gsl_vector_get(fvalues, 3) << std::endl;
#endif
    } while ((gsl_multiroot_test_residual(fvalues, tolerance) == GSL_CONTINUE) && (++niter < maxiter));
    if (niter == maxiter) {
        throw std::runtime_error("maximum iterations reached");
    }
    MArray1d costate(boost::extents[6]);
    gsl_vector *froots = gsl_multiroot_fsolver_root(solver);
    costate[0] = gsl_vector_get(froots, 0);
    costate[1] = gsl_vector_get(froots, 1);
    costate[2] = gsl_vector_get(froots, 2);
    costate[3] = gsl_vector_get(froots, 3);
    costate[4] = cop.lp[0][4];
    costate[5] = cop.lp[0][5];

    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(co_try);
    return costate;
}
