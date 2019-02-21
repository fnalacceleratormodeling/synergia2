#include "ff_rbend.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


namespace
{
    struct Rbend_adjusting_params
    {
        Rbend_adjusting_params(Lattice_element_sptr rbend, Reference_particle const & ref_part)
        : commxx_sptr(new Commxx())
        , bunch(ref_part, commxx_sptr->get_size(), 1.0e10, commxx_sptr)
        , slice(rbend)
        , ff_rbend()
        {
        }

        Commxx_sptr commxx_sptr;
        Bunch bunch;
        Lattice_element_slice slice;
        FF_rbend ff_rbend;
    };

    double rbend_propagate(double x, void * params)
    {
        struct Rbend_adjusting_params * p = (struct Rbend_adjusting_params *)params;

        Lattice_element_slice & slice = p->slice;
        Bunch & bunch = p->bunch;

        MArray2d_ref part = bunch.get_local_particles();

        part[0][0] = x;
        part[0][1] = 0;
        part[0][2] = 0;
        part[0][3] = 0;
        part[0][4] = 0;
        part[0][5] = 0;

        p->ff_rbend.apply(slice, bunch);

        //std::cout << "x0=" << x << ", x=" << part[0][0] << ", " << part[0][1] << ", " << part[0][2] << ", " << part[0][3] << ", " << part[0][4] << ", " << part[0][5] << "\n";

        return part[0][0] - x;
    }
}

double FF_rbend::adjust_rbend(Lattice_element_sptr rbend, Reference_particle const & ref_part)
{
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    gsl_function F;

    struct Rbend_adjusting_params params(rbend, ref_part);

    F.function = &rbend_propagate;
    F.params = &params;

    std::cout << F.function(-5, &params) << "\n";
    std::cout << F.function(+5, &params) << "\n";

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, -5.0, +5.0);

    int iter = 0;
    int max_iter = 100;

    double r = 0.0;
    double r0 = 0.0;

    int status;

    do
    {
        ++iter;
        status = gsl_root_fsolver_iterate(s);

        r0 = r;
        r = gsl_root_fsolver_root(s);

        double x_lo = gsl_root_fsolver_x_lower(s);
        double x_hi = gsl_root_fsolver_x_upper(s);

        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
        //status = gsl_root_test_delta(r, r0, 1e-8, 0.0);

        //std::cout << x_lo << "\t" << x_hi << "\t" << r << "\n";

    } while(status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return r;
}

FF_rbend::FF_rbend()
{

}

double FF_rbend::get_reference_cdt(double length, double strength, double angle,
                                   std::complex<double> const & phase,
                                   std::complex<double> const & term,
                                   Reference_particle &reference_particle)
{
    if (length == 0)
    {
        reference_particle.set_state_cdt(0.0);
        return 0.0;
    }
    else
    {
        double pref = reference_particle.get_momentum();
        double m = reference_particle.get_mass();

        double theta = angle / 2.0;
        double ct = cos(-theta);
        double st = sin(-theta);

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

        FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop,
                    0.0, strength, pref, m, 0.0, phase, term);

        FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

        reference_particle.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}


double FF_rbend::get_reference_cdt(double length, double angle, double edge_k_p, double * k,
                                   Reference_particle &reference_particle)
{
    if (length == 0) 
    {
        reference_particle.set_state_cdt(0.0);
        return 0.0;
    }
    else 
    {
        double pref = reference_particle.get_momentum();
        double m = reference_particle.get_mass();

        double step_length = length/steps;
        double step_strength[6] =
            { k[0]*step_length, k[1]*step_length,
              k[2]*step_length, k[3]*step_length,
              k[4]*step_length, k[5]*step_length };

        double theta = angle / 2.0;
        double ct = cos(-theta);
        double st = sin(-theta);

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        // upstream slot
        FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

        // upstream edge
        FF_algorithm::edge_unit(y, yp, edge_k_p);             // edge kick based on closed orbit
        // FF_algorithm::edge_unit(y, xp, yp, dpop, edge_k_p);    // edge kick based on particle angle

        // bend
#if 0
        FF_algorithm::yoshida<double, FF_algorithm::thin_rbend_unit<double>, 6, 3 >
            ( x, xp, y, yp, cdt, dpop,
              pref, m,
              0.0, step_length, step_strength, steps );
#endif

        // downstream edge
        FF_algorithm::edge_unit(y, yp, edge_k_p);             // edge kick based on closed orbit
        // FF_algorithm::edge_unit(y, xp, yp, dpop, -edge_k_p);   // edge kick based on particle angle

        // upstream slot
        FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

        reference_particle.set_state(x, xp, y, yp, cdt, dpop);
        return cdt;
    }
}

void FF_rbend::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF rbend on JetParticle not implemented");

#if 0
    double length = slice.get_right() - slice.get_left();
    double angle = slice.get_lattice_element().get_double_attribute("angle");
    double l     = slice.get_lattice_element().get_double_attribute("l");

    double k[6];
    k[2]  = slice.get_lattice_element().get_double_attribute("k1");
    k[3]  = 0;  // k1s
    k[4]  = slice.get_lattice_element().get_double_attribute("k2");
    k[5]  = 0;  // k2s

    typedef PropagatorTraits<JetParticle>::State_t State_t;
    typedef PropagatorTraits<JetParticle>::Component_t Component_t;

    State_t& state = jet_particle.State();

    Component_t & x(state[Chef::x]);
    Component_t & xp(state[Chef::xp]);
    Component_t & y(state[Chef::y]);
    Component_t & yp(state[Chef::yp]);
    Component_t & cdt(state[Chef::cdt]);
    Component_t & dpop(state[Chef::dpop]);

    double reference_momentum = jet_particle.ReferenceMomentum();
    double reference_brho     = jet_particle.ReferenceBRho();
    double m = jet_particle.Mass();

    k[0] = reference_brho * angle / l;
    k[1] = 0;
#endif
}

void FF_rbend::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double slice_arc = slice.get_right() - slice.get_left();

    double     a = slice.get_lattice_element().get_double_attribute("angle");
    double     l = slice.get_lattice_element().get_double_attribute("l");
    double rho_l = l / ( 2.0 * sin( a / 2.0 ) );
    double   arc = a * rho_l;

    double slice_a = a * (slice_arc / arc);
    double slice_l = l * (slice_arc / arc);

    bool ledge = slice.has_left_edge();
    bool redge = slice.has_right_edge();

    int cf = 0;
    double k_l[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // k1 and k1s
    k_l[0]  = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k_l[1]  = 0; 
    if (k_l[0] != 0.0) cf = 1;

    // k2 and k2s
    k_l[2]  = slice.get_lattice_element().get_double_attribute("k2", 0.0);
    k_l[3]  = 0;  // k2s
    if (k_l[2] != 0.0) cf = 2;

#if 0
    // k0
    k_l[0] = 2.0 * sin( a / 2.0 ) / l;
    k_l[1] = 0;

    // k1 and k1s
    k_l[2]  = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k_l[3]  = 0; 
    if (k_l[2] != 0.0) cf = 1;

    // k2 and k2s
    k_l[4]  = slice.get_lattice_element().get_double_attribute("k2", 0.0);
    k_l[5]  = 0;  // k2s
    if (k_l[4] != 0.0) cf = 2;
#endif

    // e1 and e2
    double e1 = slice.get_lattice_element().get_double_attribute("e1", 0.0);
    double e2 = slice.get_lattice_element().get_double_attribute("e2", 0.0);

    // geometries
    double theta = slice_a / 2.0;
    double ct = cos(-theta);
    double st = sin(-theta);

    double ce1 = cos(-(theta + e1));
    double se1 = sin(-(theta + e1));
    double ce2 = cos(-(theta + e2));
    double se2 = sin(-(theta + e2));

    double usAngle = theta + e1;
    double dsAngle = -(theta + e2);
    double usFaceAngle = e1;
    double dsFaceAngle = e2;

    if (!redge)
    {
        dsAngle = -theta;
        dsFaceAngle = 0;
    }

    if (!ledge)
    {
        usAngle = theta;
        usFaceAngle = 0;
    }

    double psi = - ( usFaceAngle + dsFaceAngle );
    double dphi = -psi;
    std::complex<double> phase = std::exp( std::complex<double>(0.0, psi) );
    std::complex<double> term = slice_l * std::complex<double> ( cos(dsFaceAngle), -sin(dsFaceAngle) );

    // charge, strength, and scaling
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    int charge_l = ref_l.get_charge();
    int charge_b = ref_b.get_charge();

    double m_l = ref_l.get_mass();
    double m_b = bunch.get_mass();

    double pref_l = ref_l.get_momentum();
    double pref_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double brho_l = pref_l / (PH_CNV_brho_to_p * charge_l);
    double brho_b = pref_b / (PH_CNV_brho_to_p * charge_b);

    double s_brho_l = pref_l / ref_l.get_charge();  // scaled, in GV/c
    double s_brho_b = pref_b / ref_l.get_charge();  // scaled, in GV/c

    // scaling factor (lattice to bunch)
    double scale = s_brho_l / s_brho_b;

    double k_b[6] = { 
        k_l[0] * scale, k_l[1] * scale,
        k_l[2] * scale, k_l[3] * scale,
        k_l[4] * scale, k_l[4] * scale
    };

    // magnetic field
    double strength = brho_l / rho_l;
    double eB = charge_b * strength;

    // edge strength
    double us_edge_k  =   ((charge_l > 0) ? 1.0 : -1.0 ) * strength * tan(usAngle) / brho_b;
    double ds_edge_k  = - ((charge_l > 0) ? 1.0 : -1.0 ) * strength * tan(dsAngle) / brho_b;

    double us_edge_kp =   ((charge_l > 0) ? 1.0 : -1.0 ) * strength / brho_b;
    double ds_edge_kp = - ((charge_l > 0) ? 1.0 : -1.0 ) * strength / brho_b;

    double us_edge_kx = 0.0;
    double us_edge_ky = 0.0;

    double ds_edge_kx = 0.0;
    double ds_edge_ky = 0.0;

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    if (cf==0)
    {
#if 0
        double reference_cdt = get_reference_cdt(slice_l, strength, slice_a, phase, term, ref_l);
#endif

        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        if (ledge)
        {
            // slot
            FF_algorithm::slot_unit(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, ce1, se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            us_edge_kx = us_edge_kp * (xp_l/zp_l);
            us_edge_ky = us_edge_kp * (yp_l/zp_l);

            // edge
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, us_edge_kx/scale, us_edge_ky/scale, 0);
        }

        FF_algorithm::bend_unit(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
                    0.0, strength, pref_l, m_l, 0.0, phase, term);

        if (redge)
        {
            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            ds_edge_kx = ds_edge_kp * (xp_l/zp_l);
            ds_edge_ky = ds_edge_kp * (yp_l/zp_l);

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, ds_edge_kx/scale, ds_edge_ky/scale, 0);

            // slot
            FF_algorithm::slot_unit(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, ce2, se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        double ref_cdt = cdt_l;

        // bunch particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                if (ledge)
                {
                    // upstream slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref_b, m_b);

                    // upstream edge
                    // FF_algorithm::edge_unit(y, yp, edge_k);             // edge kick based on closed orbit
                    // FF_algorithm::edge_unit(y, xp, yp, dpop, edge_k_p); // edge kick based on particle angle
                    FF_algorithm::edge_unit(y, xp, yp, us_edge_kx, us_edge_ky, 0);
                }

                // bend
                // FF_algorithm::dipole_unit(x, xp, y, yp, cdt, dpop, length, k[0]);
                FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop,
                        dphi, eB, pref_b, m_b, ref_cdt, phase, term);

                if (redge)
                {
                    // downstream edge
                    // FF_algorithm::edge_unit(y, yp, -edge_k);             // edge kick based on closed orbit
                    // FF_algorithm::edge_unit(y, xp, yp, dpop, -edge_k_p); // edge kick based on particle angle
                    FF_algorithm::edge_unit(y, xp, yp, ds_edge_kx, ds_edge_ky, 0);

                    // downstream slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref_b, m_b);
                }

                particles[part][Bunch::x]   = x;
                particles[part][Bunch::xp]  = xp;
                particles[part][Bunch::y]   = y;
                particles[part][Bunch::yp]  = yp;
                particles[part][Bunch::cdt] = cdt;
            }
        }

        // bunch spectator particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x   (s_particles[part][Bunch::x   ]);
                double xp  (s_particles[part][Bunch::xp  ]);
                double y   (s_particles[part][Bunch::y   ]);
                double yp  (s_particles[part][Bunch::yp  ]);
                double cdt (s_particles[part][Bunch::cdt ]);
                double dpop(s_particles[part][Bunch::dpop]);

                if (ledge)
                {
                    // upstream slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref_b, m_b);

                    // upstream edge
                    // FF_algorithm::edge_unit(y, yp, edge_k);             // edge kick based on closed orbit
                    // FF_algorithm::edge_unit(y, xp, yp, dpop, edge_k_p); // edge kick based on particle angle
                    FF_algorithm::edge_unit(y, xp, yp, us_edge_kx, us_edge_ky, 0);
                }

                // bend
                // FF_algorithm::dipole_unit(x, xp, y, yp, cdt, dpop, length, k[0]);
                FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop,
                        dphi, eB, pref_b, m_b, ref_cdt, phase, term);

                if (redge)
                {
                    // downstream edge
                    // FF_algorithm::edge_unit(y, yp, -edge_k);             // edge kick based on closed orbit
                    // FF_algorithm::edge_unit(y, xp, yp, dpop, -edge_k_p); // edge kick based on particle angle
                    FF_algorithm::edge_unit(y, xp, yp, ds_edge_kx, ds_edge_ky, 0);

                    // downstream slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref_b, m_b);
                }

                s_particles[part][Bunch::x]   = x;
                s_particles[part][Bunch::xp]  = xp;
                s_particles[part][Bunch::y]   = y;
                s_particles[part][Bunch::yp]  = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }
    }
    else
    {
#if 0
        // drift - kick (includes the dipole, quad, and sextupole kicks) - drift
        double step_length = slice_l/steps;
        double step_strength[6] = { k_b[0]*step_length, k_b[1]*step_length,
                                    k_b[2]*step_length, k_b[3]*step_length,
                                    k_b[4]*step_length, k_b[5]*step_length };

        double reference_cdt = get_reference_cdt(slice_l, slice_a, edge_k/scale, k_l, ref_l);
        // double reference_cdt = get_reference_cdt(length, angle, edge_k_p/scale, k, ref_l);
        double step_reference_cdt = reference_cdt/steps;

        // with combined high order function, use yoshida approximation
        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part)
        {
            double x   (particles[part][Bunch::x   ]);
            double xp  (particles[part][Bunch::xp  ]);
            double y   (particles[part][Bunch::y   ]);
            double yp  (particles[part][Bunch::yp  ]);
            double cdt (particles[part][Bunch::cdt ]);
            double dpop(particles[part][Bunch::dpop]);

            // upstream slot
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref_b, m);

            // upstream edge
            FF_algorithm::edge_unit(y, yp, edge_k);             // edge kick based on closed orbit
            //FF_algorithm::edge_unit(y, xp, yp, dpop, edge_k_p);    // edge kick based on particle angle

            // bend
            FF_algorithm::yoshida<double, FF_algorithm::thin_rbend_unit<double>, 6, 3 >
                ( x, xp, y, yp, cdt, dpop,
                  pref_b, m,
                  step_reference_cdt,
                  step_length, step_strength, steps );

            // downstream edge
            FF_algorithm::edge_unit(y, yp, edge_k);             // edge kick based on closed orbit
            // FF_algorithm::edge_unit(y, xp, yp, dpop, -edge_k_p);   // edge kick based on particle angle

            // downstream slot
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref_b, m);

            particles[part][Bunch::x]  = x;
            particles[part][Bunch::xp] = xp;
            particles[part][Bunch::y]  = y;
            particles[part][Bunch::yp] = yp;
            particles[part][Bunch::cdt] = cdt;
        }
#endif

        // bend - kick (quad and sextupole) - bend
        std::complex<double> phase_e1 = FF_algorithm::bend_edge_phase(e1);
        std::complex<double> phase_e2 = FF_algorithm::bend_edge_phase(e2);

        double step_angle  = slice_a / steps;
        double step_length = slice_l / steps;

        double step_kl_l[6] = 
                { k_l[0] * step_length, k_l[1] * step_length,
                  k_l[2] * step_length, k_l[3] * step_length,
                  k_l[4] * step_length, k_l[5] * step_length };

        double step_kl_b[6] = 
                { k_b[0] * step_length, k_b[1] * step_length,
                  k_b[2] * step_length, k_b[3] * step_length,
                  k_b[4] * step_length, k_b[5] * step_length };

#if 0
        double ref_cdt = get_reference_cdt(length, angle, strength, ledge, redge,
                e1, e2, us_edge_k_p/scale, ds_edge_k_p/scale, phase_e1, phase_e2, kl, ref_l);
#endif

        // propagate the reference particle, and set the edge kick strength 
        // from the reference particle
        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        if (ledge)
        {
            // slot
            FF_algorithm::slot_unit(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, ce1, se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            us_edge_kx = us_edge_kp * (xp_l/zp_l);
            us_edge_ky = us_edge_kp * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, us_edge_k/scale);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, us_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, us_edge_kx/scale, us_edge_ky/scale, 0);

            // bend edge (thin)
            FF_algorithm::bend_edge(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, e1, phase_e1, strength, pref_l, m_l);
        }

        FF_algorithm::bend_yoshida6< double, 
                                     FF_algorithm::rbend_thin_cf_kick<double>, 
                                     FF_algorithm::rbend_unit_phase,
                                     FF_algorithm::rbend_unit_term,
                                     FF_algorithm::rbend_dphi,
                                     2 >
            ( x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
              pref_l, m_l, 0.0 /* step ref_cdt */,
              step_length, step_kl_l,
              rho_l, strength, steps );

        if (redge)
        {
            // bend edge (thin)
            FF_algorithm::bend_edge(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, e2, phase_e2, strength, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            ds_edge_kx = ds_edge_kp * (xp_l/zp_l);
            ds_edge_ky = ds_edge_kp * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, ds_edge_k);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, ds_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, ds_edge_kx/scale, ds_edge_ky/scale, 0);

            // slot
            FF_algorithm::slot_unit(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, ce2, se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        double ref_cdt = cdt_l;
        double step_ref_cdt = ref_cdt / steps;

        // with combined function
        // bunch particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                if (ledge)
                {
                    // slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref_b, m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, us_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                    FF_algorithm::edge_unit(y, xp, yp, us_edge_kx, us_edge_ky, 0);

                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e1, phase_e1, strength, pref_b, m_b);
                }

                // bend body
                FF_algorithm::bend_yoshida6< double, 
                                             FF_algorithm::rbend_thin_cf_kick<double>, 
                                             FF_algorithm::rbend_unit_phase,
                                             FF_algorithm::rbend_unit_term,
                                             FF_algorithm::rbend_dphi,
                                             2 >
                    ( x, xp, y, yp, cdt, dpop,
                      pref_b, m_b, step_ref_cdt,
                      step_length, step_kl_b,
                      rho_l, strength, steps );

                if (redge)
                {
                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e2, phase_e2, strength, pref_b, m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    FF_algorithm::edge_unit(y, xp, yp, ds_edge_kx, ds_edge_ky, 0);

                    // slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref_b, m_b);
                }

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }
        }

        // bunch spectator particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x   (s_particles[part][Bunch::x   ]);
                double xp  (s_particles[part][Bunch::xp  ]);
                double y   (s_particles[part][Bunch::y   ]);
                double yp  (s_particles[part][Bunch::yp  ]);
                double cdt (s_particles[part][Bunch::cdt ]);
                double dpop(s_particles[part][Bunch::dpop]);

                if (ledge)
                {
                    // slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref_b, m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, us_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                    FF_algorithm::edge_unit(y, xp, yp, us_edge_kx, us_edge_ky, 0);

                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e1, phase_e1, strength, pref_b, m_b);
                }

                // bend body
                FF_algorithm::bend_yoshida6< double, 
                                             FF_algorithm::rbend_thin_cf_kick<double>, 
                                             FF_algorithm::rbend_unit_phase,
                                             FF_algorithm::rbend_unit_term,
                                             FF_algorithm::rbend_dphi,
                                             2 >
                    ( x, xp, y, yp, cdt, dpop,
                      pref_b, m_b, step_ref_cdt,
                      step_length, step_kl_b,
                      rho_l, strength, steps );

                if (redge)
                {
                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e2, phase_e2, strength, pref_b, m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    FF_algorithm::edge_unit(y, xp, yp, ds_edge_kx, ds_edge_ky, 0);

                    // slot
                    FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref_b, m_b);
                }

                s_particles[part][Bunch::x]  = x;
                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::y]  = y;
                s_particles[part][Bunch::yp] = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }
    }

    bunch.get_reference_particle().increment_trajectory(slice_l);
}

template<class Archive>
    void
    FF_rbend::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_rbend::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_rbend::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_rbend::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_rbend::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_rbend::~FF_rbend()
{

}

