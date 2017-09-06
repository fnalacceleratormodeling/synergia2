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
#include <iomanip>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


#define DEBUG 0
#define TRYDEBUG 0

struct Closed_orbit_params
{
    const double dpp;
    Lattice_sptr working_lattice_sptr;
    Stepper_sptr stepper_sptr;
    Closed_orbit_params(const double dpp, Lattice_sptr lattice_sptr) :
        dpp(dpp)
        , working_lattice_sptr(new Lattice(*lattice_sptr))
    {
#if DEBUG
        std::cout << "egs: Closed_orbit_params constructor, dpp" << dpp << std::endl;
        std::cout << "egs: lattice length: " << working_lattice_sptr->get_length() << std::endl;
#endif
        // turn off any RF cavities because they screw up the closed orbit calcation
        for (Lattice_elements::const_iterator lep=working_lattice_sptr->get_elements().begin();
             lep!=working_lattice_sptr->get_elements().end(); ++lep) {
            if ((*lep)->get_type() == "rfcavity") {
                (*lep)->set_double_attribute("volt", 0.0);
            }
        }
        stepper_sptr = Stepper_sptr(new Independent_stepper(working_lattice_sptr, 1, 1));
    }
};

int
propagate_co_try(const gsl_vector *co_try, void *params, gsl_vector *co_results)
{
    Closed_orbit_params *copp = static_cast<Closed_orbit_params *>(params);

    Commxx_sptr commxx_sptr(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(copp->working_lattice_sptr->get_reference_particle(),
                                   commxx_sptr->get_size(), 1.0e10, commxx_sptr));
    MArray2d_ref lp(bunch_sptr->get_local_particles());
    Bunch_simulator bunch_simulator(bunch_sptr);

    Propagator propagator(copp->stepper_sptr);
    if (bunch_sptr->get_total_num() == 0) {
        throw std::runtime_error("Lattice will not transport particles.  Is it unstable?");
    }

    // set phase space coordinates in bunch
    double orbit_start[4];
    orbit_start[0] = lp[0][0] = gsl_vector_get(co_try, 0);
    orbit_start[1] = lp[0][1] = gsl_vector_get(co_try, 1);
    orbit_start[2] = lp[0][2] = gsl_vector_get(co_try, 2);
    orbit_start[3] = lp[0][3] = gsl_vector_get(co_try, 3);
    lp[0][4] = 0.0;
    lp[0][5] = copp->dpp;
    
#if TRYDEBUG
    std::cout << "egs: propagate_co_try: initial: " << std::setprecision(15) << lp[0][0]<<", "<<lp[0][1]<<", "<<lp[0][2]<<", "<<lp[0][3]<<std::endl;
#endif
    propagator.propagate(bunch_simulator, 1, 1, DEBUG);
#if TRYDEBUG
    std::cout << "egs: propagate_co_try: final: " << std::setprecision(15) << lp[0][0]<<", "<<lp[0][1]<<", "<<lp[0][2]<<", "<<lp[0][3]<<std::endl;
#endif
    gsl_vector_set(co_results, 0, lp[0][0] - orbit_start[0]);
    gsl_vector_set(co_results, 1, lp[0][1] - orbit_start[1]);
    gsl_vector_set(co_results, 2, lp[0][2] - orbit_start[2]);
    gsl_vector_set(co_results, 3, lp[0][3] - orbit_start[3]);
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
    //const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid;
    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrids;
    //const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_dnewton;
    gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(T, ndim);
#if DEBUG
    std::cout << "egs: after gsl_multiroot_fsolver_alloc" << std::endl;
#endif
    // co_try are the coordinates of the closed orbit
    gsl_vector *co_try = gsl_vector_alloc(ndim);
    // initialize the closed orbit
    for (int i=0; i<ndim; ++i) {
        gsl_vector_set(co_try, i, 0.0);
    }
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
            throw std::runtime_error("Closed orbit solver unable to converge.  Is the tolerance too tight?");
            break;
        case GSL_EBADFUNC:
#if DEBUG
            std::cout << "egs: EBADFUNC" << std::endl;
#endif
            throw std::runtime_error("Closed orbit solver failed to evaluate solution");
            break;
        default:
#if DEBUG
            std::cout << "egs: success" << std::endl;
#endif
            break;
        }
        fvalues = gsl_multiroot_fsolver_f(solver);
#if DEBUG
        std::cout << "egs: residuals at current iteration: " << gsl_vector_get(solver->f, 0) << ", "
            << gsl_vector_get(solver->f, 1) << ", " << gsl_vector_get(solver->f, 2) <<
            ", " << gsl_vector_get(solver->f, 3) << std::endl;
        gsl_vector *froots = gsl_multiroot_fsolver_root(solver);
        std::cout << "egs: roots at current iteration: " << gsl_vector_get(solver->x, 0) << ", "
            << gsl_vector_get(solver->x, 1) << ", " << gsl_vector_get(solver->x, 2) <<
            ", " << gsl_vector_get(solver->x, 3) << std::endl;
#endif
    } while ((gsl_multiroot_test_residual(solver->f, tolerance) == GSL_CONTINUE) && (++niter < maxiter));
    if (niter == maxiter) {
        std::stringstream sstr;
        sstr << "Could not locate closed orbit after " << maxiter << " iterations";
        throw std::runtime_error( sstr.str() );
    }
    MArray1d costate(boost::extents[6]);
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
