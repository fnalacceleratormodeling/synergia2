
#ifdef __CUDA_ARCH__

// no implementations for CUDA arch

#else


#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"

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
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, Commxx().size(), 1e09);

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
        bunch.checkout_particles();

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

        gsl_multiroot_fsolver_f(solver);

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

#include "synergia/foundation/trigon.h"

std::array<double, 2>
filter_transverse_tunes(double const* jac);
 
// [tune_h, tune_v, c_delta_t]
std::array<double, 3>
Lattice_simulator::calculate_tune_and_cdt(Lattice const& lattice, double dpp)
{
    // trigon bunch
    using trigon_t = Trigon<double, 2, 6>;

    // get the reference particle
    auto const& ref = lattice.get_reference_particle();

    // closed orbit
    auto probe = Lattice_simulator::calculate_closed_orbit(lattice, dpp);

    // comm world
    Commxx comm;

    bunch_t<trigon_t> tb(ref, comm.size(), comm);
    bunch_t<double>   pb(ref, comm.size(), 1e9, comm);

    // design reference particle from the closed orbit
    auto ref_l = ref;
    ref_l.set_state(probe);
    tb.set_design_reference_particle(ref_l);
    pb.set_design_reference_particle(ref_l);

    auto tparts = tb.get_host_particles();
    auto pparts = pb.get_host_particles();

    // init value
    for(int i=0; i<6; ++i) 
    {
        tparts(0, i).set(probe[i], i);
        pparts(0, i) = probe[i];
    }

    // check in
    tb.checkin_particles();
    pb.checkin_particles();

    // init c_delta_t
    double c_delta_t = 0.0;

    // propagate trigon
    for(auto & ele : lattice.get_elements())
    {
        if (ele.get_type() == element_type::rfcavity)
        {
            Lattice_element dup = ele;
            dup.set_double_attribute("volt", 0.0);

            FF_element::apply(dup, tb);
            FF_element::apply(dup, pb);
        }
        else
        {
            FF_element::apply(ele, tb);
            FF_element::apply(ele, pb);
        }

        // cdt from reference particle
        c_delta_t += pb.get_design_reference_particle()
                       .get_state()[Bunch::cdt];
    }

    // checkout particles
    tb.checkout_particles();
    pb.checkout_particles();

    // cdt from actual particle
    c_delta_t += pparts(0, 4);

    // one-turn-map
    auto kjac = tb.get_jacobian(0);

    //auto jac = karray_to_matrix(kjac);
    //auto nus = filter_transverse_tunes(jac);
    auto nus = filter_transverse_tunes(kjac.data());

#if 0
    std::cout << "nus0 = " << nus[0] << "\n";
    std::cout << "nus1 = " << nus[1] << "\n";
    std::cout << "c_delta_t = " << c_delta_t << "\n";
#endif

    double tune_h = nus[0];
    double tune_v = nus[1];

    return {tune_h, tune_v, c_delta_t};
}


chromaticities_t
Lattice_simulator::get_chromaticities(Lattice const& lattice, double dpp)
{
    chromaticities_t chroms;

    auto ref = lattice.get_reference_particle();
    double gamma = ref.get_gamma();

    // tune = [tune_h, tune_v, cdt]
    auto tune_0  = calculate_tune_and_cdt(lattice, 0.0);
    auto tune_p  = calculate_tune_and_cdt(lattice, dpp);
    auto tune_m  = calculate_tune_and_cdt(lattice, -dpp);
    auto tune_pp = calculate_tune_and_cdt(lattice, 2.0*dpp);
    auto tune_mm = calculate_tune_and_cdt(lattice, -2.0*dpp);

    // five point stencil:
    // given function f(x) = a_0 + a_1*x + a_2*x**2 + a_3*x**3 + a_4*x**4
    // choose offset spacing d, calculate values
    // y++ = f(2*d)
    // y+ = f(d)
    // y0 = f(0)
    // y- = f(-d)
    // y-- = f(-2*d)
    // then you can easily calculate:
    // y0 = a_0
    // y+ - y- = 2*a_1*d + 2*a_3*d**3
    // y++ - y-- = 4*a_1*d + 16*a_3*d**3
    // y+ + y- - 2*y0 = 2*a_2*d**2 + 2*a_4*d**4
    // y++ + y-- - 2*y0 = 8*a_2*d**2 + 32*a_4*d**4
    // yielding expressions:
    // 8*(y+ - y-) - (y++ - y--) = 12*a_1*d
    // (y++ - y--) - 2*(y+ - y-) = 12*a_3*d**3
    // 16*(y+ + y- - 2*y_0) - (y++ + y-- - 2*y_0) = 24*a_2*d**2
    // (y++ + y-- - 2*y_0) - 4*(y+ + y- - 2*y_0) = 24*a_4*d**4

    double tune_h0 = tune_0[0];
    double tune_v0 = tune_0[1];
    double cT0 = tune_0[2];

    double tune_h_plus = tune_p[0];
    double tune_v_plus = tune_p[1];
    double c_delta_t_plus = tune_p[2];

    double tune_h_minus = tune_m[0];
    double tune_v_minus = tune_m[1];
    double c_delta_t_minus = tune_m[2];

    double tune_h_plusplus = tune_pp[0];
    double tune_v_plusplus = tune_pp[1];
    double c_delta_t_plusplus = tune_pp[2];

    double tune_h_minusminus = tune_mm[0];
    double tune_v_minusminus = tune_mm[1];
    double c_delta_t_minusminus = tune_mm[2];

#if 0
    std::cout << "tune_0: " << tune_0[0] << ", " << tune_0[1] << ", " << tune_0[2] << "\n";
    std::cout << "tune_p: " << tune_p[0] << ", " << tune_p[1] << ", " << tune_p[2] << "\n";
    std::cout << "tune_m: " << tune_m[0] << ", " << tune_m[1] << ", " << tune_m[2] << "\n";
    std::cout << "tune_pp: " << tune_pp[0] << ", " << tune_pp[1] << ", " << tune_pp[2] << "\n";
    std::cout << "tune_mm: " << tune_mm[0] << ", " << tune_mm[1] << ", " << tune_mm[2] << "\n";
#endif

    double a_h_chrom, b_h_chrom;
    a_h_chrom = 0.5*(tune_h_plus-tune_h_minus)/dpp;
    b_h_chrom = 0.25*(tune_h_plusplus-tune_h_minusminus)/dpp;
    double horizontal_chromaticity_alt=(4*a_h_chrom - b_h_chrom)/3.0;

    double a_v_chrom, b_v_chrom;
    a_v_chrom = 0.5*(tune_v_plus-tune_v_minus)/dpp;
    b_v_chrom = 0.25*(tune_v_plusplus-tune_v_minusminus)/dpp;
    double vertical_chromaticity_alt=(4*a_v_chrom - b_v_chrom)/3.0;

    double a_slip, b_slip;
    a_slip = 0.5*(c_delta_t_plus-c_delta_t_minus)/cT0 / dpp;
    b_slip = 0.25*(c_delta_t_plusplus-c_delta_t_minusminus)/cT0 / dpp;
    double slip_factor_alt =(4*a_slip-b_slip)/3.0;
    double cdtp_m_cdtm = c_delta_t_plus - c_delta_t_minus;
    double cdtp_p_cdtm_m2cdt0 = c_delta_t_plus + c_delta_t_minus - 2*cT0;
    double cdtpp_m_cdtmm = c_delta_t_plusplus - c_delta_t_minusminus;
    double cdtpp_p_cdtmm_m2cdt0 = c_delta_t_plusplus + c_delta_t_minusminus - 2*cT0;

    double qxp_m_qxm = tune_h_plus - tune_h_minus;
    double qxp_p_qxm_m2qx0 = tune_h_plus + tune_h_minus - 2*tune_h0;
    double qxpp_m_qxmm = tune_h_plusplus - tune_h_minusminus;
    double qxpp_p_qxmm_m2qx0 = tune_h_plusplus + tune_h_minusminus - 2*tune_h0;

    double qyp_m_qym = tune_v_plus - tune_v_minus;
    double qyp_p_qym_m2qy0 = tune_v_plus + tune_v_minus - 2*tune_v0;
    double qypp_m_qymm = tune_v_plusplus - tune_v_minusminus;
    double qypp_p_qymm_m2qy0 = tune_v_plusplus + tune_v_minusminus - 2*tune_v0;

    chroms.slip_factor 
        = (8.0*cdtp_m_cdtm - cdtpp_m_cdtmm)/(12.0*dpp*cT0);

    chroms.slip_factor_prime 
        = 2.0*(16*cdtp_p_cdtm_m2cdt0 - cdtpp_p_cdtmm_m2cdt0) /(24.0*dpp*dpp*cT0);

    chroms.momentum_compaction 
        = chroms.slip_factor + 1. / gamma / gamma;

    // d^2 cdt/d dpop^2 = 2! *a_2

    // d^3 cdt/ d dpop^3 = 3! * a_3 but I'm not planning on using it
    // double d3cdt_d_dpop3 = 6.0*(cdtpp_m_cdtmm - 2.0*cdtp_m_cdtm)/(12.0*dpp*dpp*dpp*cT0);

    // d^4 cdt / d dpop^4 = 4! * a_4 but I'm not planning on using it
    // double d4cdt_d_dpop4 = 24.0*(cdtpp_p_cdtmm_m2cdt0 - 4.0*cdtp_p_cdtm_m2cdt0)/(24.0*dpp*dpp*dpp*dpp*cT0);

    chroms.horizontal_chromaticity 
        = (8.0*qxp_m_qxm - qxpp_m_qxmm)/(12.0*dpp);

    chroms.vertical_chromaticity 
        = (8.0*qyp_m_qym - qypp_m_qymm)/(12.0*dpp);

    chroms.horizontal_chromaticity_prime 
        = 2.0*(16*qxp_p_qxm_m2qx0 - qxpp_p_qxmm_m2qx0) /(24.0*dpp*dpp);

    chroms.vertical_chromaticity_prime 
        = 2.0*(16*qyp_p_qym_m2qy0 - qypp_p_qymm_m2qy0) /(24.0*dpp*dpp);


    return chroms;
}

karray2d_row
Lattice_simulator::get_linear_one_turn_map(Lattice const& lattice)
{
    // 2nd order one-turn-map is sufficient for the jacobian
    auto map = get_one_turn_map<2>(lattice);
    return map.jacobian();
}

std::array<double, 3>
Lattice_simulator::map_to_twiss(karray2d_row map)
{
    // map must be 2x2
    if (map.extent(0)!=2 || map.extent(1)!=2)
        throw std::runtime_error("map_to_twiss: wrong dimensions");

    // [alpha, beta, psi]
    std::array<double, 3> ret;

    double cosmu = 0.5 * (map(0, 0) + map(1, 1));
    double asinmu = 0.5 * (map(0, 0) - map(1, 1));

    if (abs(cosmu) > 1.0)
        throw std::runtime_error("map_to_twiss: map is unstable");

    double mu = acos(cosmu);

    // beta is positive
    if (map(0,1) < 0.0) mu = 2.0*mconstants::pi - mu;

    ret[0] = asinmu / sin(mu);
    ret[1] = map(0,1) / sin(mu);
    ret[2] = mu / (2.0 * mconstants::pi);

    return ret;
}

double
Lattice_simulator::get_bucket_length(Lattice const& lattice)
{
    double freq = 0.0;
    double freq2 = 0.0;

    double harmon = 0.0;
    double harmon2 = 0.0;

    bool iswf = false;
    bool iswh = false;

    double eps = 1e-6;

    for(auto & ele : lattice.get_elements())
    {
        if (ele.get_type() == element_type::rfcavity)
        {
            if (ele.has_double_attribute("harmon"))
            {
                harmon = ele.get_double_attribute("harmon");

                if (iswh && (abs(harmon-harmon2)>eps))
                    throw std::runtime_error("get_bucket_length:" 
                            " rf elements with different harmonic" 
                            " number found!");

                harmon2 = harmon;
                iswh = true;
            }

            if (ele.has_double_attribute("freq"))
            {
                freq = ele.get_double_attribute("freq");

                if (iswf && (abs(freq-freq2)>eps))
                    throw std::runtime_error("get_bucket_length:" 
                            " rf elements with different frequency" 
                            " found!");

                freq2 = freq;
                iswf = true;
            }
        }
    }

    // use harmonic number
    if (iswh) return lattice.get_length()/harmon;

    // use frequency
    if (iswf)
    {
        double beta = lattice.get_reference_particle().get_beta();
        return pconstants::c * beta / freq;
    }

    // or return 0
    return 0.0;
}


#endif // __CUDA_ARCH__

