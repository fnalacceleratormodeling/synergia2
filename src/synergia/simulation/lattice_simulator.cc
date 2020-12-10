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
    auto sim = Bunch_simulator::create_single_bunch_simulator(ref, 1, 1e09);

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

#include "synergia/utils/multi_array_conversions.h"
#include "synergia/foundation/trigon.h"

namespace
{
    std::array<double, 2>
    filter_transverse_tunes(MatrixD const& jac)
    {
        std::array<double, 2> nu;

        const int ix  = 0;
        const int ipx = 1;
        const int iy  = 2;
        const int ipy = 3;

        if ( jac(iy , ix ) || jac(ix , iy ) ||
             jac(ix , ipy) || jac(iy , ipx) ||
             jac(ipy, ix ) || jac(ipx, iy ) ||
             jac(ipy, ipx) || jac(ipx, ipy) )
        {
            auto lambda = jac.eigenvalues();

            for(int i=0; i<6; ++i) 
            {
                if( fabs( abs(lambda(i)) - 1.0 ) > 1.0e-4 ) 
                {
                    std::stringstream ss;
                    ss << "filterTransverseTunes: "
                       << "the lattice is nearly unstable. "
                       << "lambda( " << i << " ) has magnitude = "
                       << abs(lambda(i));
                    throw std::runtime_error(ss.str());
                }
            }

            if( (abs(lambda(0) - std::conj(lambda(1))) > 1.0e-4)  ||
                (abs(lambda(3) - std::conj(lambda(4))) > 1.0e-4) ) 
            {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "conjugacy condition has been vilated. "
                   << "The lattice may be linearly unstable. "
                   << "Eigenvalues = " << lambda;
                throw std::runtime_error(ss.str());
            }
           
            double csH = lambda(0).real();
            double csV = lambda(3).real();

            if( fabs( csH - csV ) < 1.0e-4 ) 
            {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "\"Horizontal\" and \"vertical\" tunes "
                   << "are too near each other for reasonable results. "
                   << "The calculation is meaningless.";
                throw std::runtime_error(ss.str());
            }

            double  dcos, cos2phi, sin2phi, tanphi;

            MatrixD U( 2, 2 ), S( 2, 2 );

            U << 1.0, 0.0,  0.0, 1.0;
            S << 0.0, 1.0, -1.0, 0.0;

            MatrixD M( 2, 2 ), N( 2, 2 );
            MatrixD m( 2, 2 ), n( 2, 2 );

            M << jac(ix,ix), jac(ix,ipx), jac(ipx,ix), jac(ipx,ipx);
            N << jac(iy,iy), jac(iy,ipy), jac(ipy,iy), jac(ipy,ipy);
            m << jac(iy,ix), jac(iy,ipx), jac(ipy,ix), jac(ipy,ipx);
            n << jac(ix,iy), jac(ix,ipy), jac(ipx,iy), jac(ipx,ipy);

#if 0
            std::cout << "M = " << M << "\n";
            std::cout << "N = " << N << "\n";
            std::cout << "m = " << m << "\n";
            std::cout << "n = " << n << "\n";
#endif

            dcos    = csH - csV;
            cos2phi = ( M - N ).trace() / ( 2.0 *( dcos ) );

            if( fabs(cos2phi - 1.0) < 1.0e-4 ) 
                cos2phi =   1.0;  // ??? Rather coarse,

            if( fabs(cos2phi + 1.0) < 1.0e-4 ) 
                cos2phi = - 1.0;  // ??? isn't it?
            
            if( fabs(cos2phi) > 1.0 ) 
            {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "cos(2 phi) = " << std::setprecision(10) << cos2phi 
                   << "; has magnitude larger than one. "
                   << "Cannot continue calculation. ";
                throw std::runtime_error(ss.str());
            }
            
            if( cos2phi < 0.0 ) 
            {
                sin2phi = csH;  // Variable used as dummy register.
                csH     = csV;
                csV     = sin2phi;
                dcos    = -dcos;
                cos2phi = -cos2phi;
            }

            sin2phi = sqrt( 1.0 - cos2phi*cos2phi );
            tanphi  = sin2phi / ( 1.0 + cos2phi );

#if 0
            std::cout << "sin2phi = " << sin2phi << "\n";
            std::cout << "tanphi = " << tanphi << "\n";
#endif
            
            MatrixD D( 2, 2 ), A( 2, 2 ), B( 2, 2 );

            if( fabs(sin2phi) > 1.0e-8 ) 
            {
                D = -(m + S*n.transpose()*S.transpose()) * 
                     (1.0 / (dcos*sin2phi)); 
            }
            else 
            {
                D << 1.0, 0.0, 0.0, 1.0;
            }
            
            if( fabs(D.determinant() - 1.0) > 1.0e-4 ) 
            {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "The matrix D is non-symplectic. "
                   << "|D| = " << D.determinant();
                throw std::runtime_error(ss.str());
            }
            
            // ...... Edwards-Teng sign convention.
            if( D.trace() < 0.0 ) 
            {
                D = -D;
                sin2phi = -sin2phi;
                tanphi  = -tanphi;
            }

            A = M - D.inverse()*m*tanphi;
            B = N + D*n*tanphi;

#if 0
            std::cout << "A = " << A << "\n";
            std::cout << "B = " << B << "\n";
#endif
           
            // ......  First the "horizontal" ......
            MatrixD JH = A - csH*U;
            double snH = (JH(0,1)>0.0) ?  sqrt(1.0 - csH*csH) 
                                       : -sqrt(1.0 - csH*csH);

            // .......... A little test to keep everyone honest .....
            if( JH(0,0) )
            {
                if( fabs((JH(0,0) + JH(1,1)) / (JH(0,0) - JH(1,1))) > 1.0e-4 ) 
                {
                    std::cout
                       << "WARNING -- filterTransverseTunes: "
                       << "\"Horizontal\" matrix does not "
                       << "pass symplecticity test. "
                       << "JH( 0, 0 ) = " << JH( 0, 0 ) << ", "
                       << "JH( 1, 1 ) = " << JH( 1, 1 ) << ". "
                       << "The ratio is " 
                       << fabs((JH(0,0) + JH(1,1)) / (JH(0,0) - JH(1,1)))
                       << "\n";
                }
            }
           
           
            // ......  Then  the "vertical" ......
            MatrixD JV = B - csV*U;
            double snV = (JV(0,1)>0.0) ?  sqrt(1.0 - csV*csV) 
                                       : -sqrt(1.0 - csV*csV);
           
            // .......... A little test to keep everyone honest .....
            if( JV(0,0) )
            {
                if( fabs((JV(0,0) + JV(1,1)) / (JV(0,0) - JV(1,1))) > 1.0e-4 ) 
                {
                    std::cout
                       << "WARNING -- filterTransverseTunes: "
                       << "\"Vertical\" matrix does not "
                       << "pass symplecticity test. "
                       << "JV( 0, 0 ) = " << JV( 0, 0 ) << ", "
                       << "JV( 1, 1 ) = " << JV( 1, 1 ) << ". "
                       << "The ratio is " 
                       << fabs((JV(0,0) + JV(1,1)) / (JV(0,0) - JV(1,1)))
                       << "\n";
                }
            }

            const double M_TWOPI = mconstants::pi * 2;
           
            double theta = atan2( snH, csH );
            if( theta < 0.0 ) theta += M_TWOPI;
            nu[0] = theta / M_TWOPI;

            theta = atan2( snV, csV );
            if( theta < 0.0 )  theta += M_TWOPI;
            nu[1] = theta / M_TWOPI;
        }
        else
        {
            double sn, cs;

            // Uncoupled calculation .....
            // (Lifted from LattFuncSage) ...
            // ... first horizontal
            cs = (jac(ix, ix) + jac(ipx, ipx)) / 2.0;

            if( fabs(cs) <= 1.0 ) 
            { 
                if( jac(ix, ipx) > 0.0 )  sn =   sqrt(1.0 - cs*cs);
                else                      sn = - sqrt(1.0 - cs*cs);
            }
            else 
            {
                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "cos( psi_H ) = " << cs << ". "
                   << "Cannot continue with calculation.";
                throw std::runtime_error(ss.str());
            }

            const double M_TWOPI = mconstants::pi * 2;

            double theta = atan2(sn, cs);
            if( theta < 0.0 )  theta += M_TWOPI;
            nu[0] = theta / M_TWOPI;
     

            // ... then vertical.
            cs = (jac(iy, iy) + jac(ipy, ipy)) / 2.0;

            if( fabs(cs) <= 1.0 ) 
            {
                if( jac(iy, ipy) > 0.0 )  sn =   sqrt(1.0 - cs*cs);
                else                      sn = - sqrt(1.0 - cs*cs);
            }
            else 
            {

                std::stringstream ss;
                ss << "filterTransverseTunes: "
                   << "cos( psi_V ) = " << cs << ". "
                   << "Cannot continue with calculation.";
                throw std::runtime_error(ss.str());
            }
          
            theta = atan2(sn, cs);
            if( theta < 0.0 )   theta += M_TWOPI;
            nu[1] = theta / M_TWOPI;
        }

        return nu;
    }
}

// [tune_h, tune_v, c_delta_t]
std::array<double, 3>
Lattice_simulator::calculate_tune_and_cdt(Lattice const& lattice, double dpp)
{
    // get the reference particle
    auto const& ref = lattice.get_reference_particle();

    // closed orbit
    auto probe = Lattice_simulator::calculate_closed_orbit(lattice, dpp);

#if 0
    std::cout << "closed orbit: ";
    for(int i=0; i<6; ++i)
        std::cout << probe[i] << ", ";
    std::cout << "\n";
#endif

    // trigon bunch
    using trigon_t = Trigon<double, 1, 6>;

    bunch_t<trigon_t> tb(ref);
    bunch_t<double>   pb(ref, 1, 1e9);

    auto tparts = tb.get_host_particles();
    auto pparts = pb.get_host_particles();

    // init value
    for(int i=0; i<6; ++i) 
    {
        tparts(0, i) = trigon_t(probe[i], i);
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
    auto jac = karray_to_matrix(kjac);

#if 0
    for(int i=0; i<6; ++i)
    {
        for(int j=0; j<6; ++j) std::cout << jac(i, j) << ", ";
        std::cout << "\n";
    }
#endif

    auto nus = filter_transverse_tunes(jac);

#if 0
    std::cout << "nus0 = " << nus[0] << "\n";
    std::cout << "nus1 = " << nus[1] << "\n";
    std::cout << "c_delta_t = " << c_delta_t << "\n";
#endif

    double tune_h = nus[0];
    double tune_v = nus[1];

    return {tune_h, tune_v, c_delta_t};
}



