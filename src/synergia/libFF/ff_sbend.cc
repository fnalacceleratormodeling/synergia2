#include "ff_sbend.h"

#include "synergia/foundation/physical_constants.h"

// p [Gev/c] = -- * B*rho [ Tesla meters ]
#define PH_CNV_brho_to_p   (1.0e-9 * pconstants::c)


#if 0
double FF_sbend::get_reference_cdt(double length, double strength, double angle,
                                   bool ledge, bool redge,
                                   double e1, double e2, 
                                   double us_edge_k_p, double ds_edge_k_p,
                                   double dphi,
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

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double ce1 = cos(-e1);
        double se1 = sin(-e1);
        double ce2 = cos(-e2);
        double se2 = sin(-e2);

        if (ledge)
        {
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref, m);
            FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
        }

        FF_algorithm::bend_complete(x, xp, y, yp, cdt, dpop,
                   dphi, strength, pref, m, 0.0/*ref cdt*/, phase, term);

        if (redge)
        {
            FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref, m);
        }

        reference_particle.set_state(x, xp, y, yp, cdt, dpop);
        return cdt;
    }
}


double FF_sbend::get_reference_cdt(double length, double angle, double strength,
                                   bool ledge, bool redge,
                                   double e1, double e2, 
                                   double us_edge_k_p, double ds_edge_k_p,
                                   std::complex<double> phase_e1, 
                                   std::complex<double> phase_e2, 
                                   double * kl,
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

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double ce1 = cos(-e1);
        double se1 = sin(-e1);
        double ce2 = cos(-e2);
        double se2 = sin(-e2);

        double step_length = length / steps;
        double step_angle  = angle / steps;
        double r0 = length / angle;

        double step_strength[6] = { kl[0] * step_length, kl[1] * step_length,
                                    kl[2] * step_length, kl[3] * step_length,
                                    kl[4] * step_length, kl[5] * step_length };

        if (ledge)
        {
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref, m);

            //FF_algorithm::edge_unit(y, yp, us_edge_k);
            FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);

            FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e1, phase_e1, strength, pref, m);
        }

        // bend
#if 0
        FF_algorithm::bend_yoshida6<double, FF_algorithm::thin_cf_kick_2<double>, 2>
            ( x, xp, y, yp, cdt, dpop,
              pref, m, 0.0 /* step ref_cdt */,
              step_angle, step_strength,
              r0, strength, steps );
#endif

        if (redge)
        {
            FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e2, phase_e2, strength, pref, m);

            //FF_algorithm::edge_unit(y, yp, ds_edge_k);
            FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);

            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref, m);
        }

        reference_particle.set_state(x, xp, y, yp, cdt, dpop);
        return cdt;
    }
}
#endif

namespace
{
    struct SbendParams
    {
        bool ledge;
        bool redge;

        double length;
        double angle;
        double r0;

        int    steps;
        double step_kl[6];

        double dphi;
        double strength;
        double pref_b;
        double m_b;
        double ref_cdt;

        double e1;
        double e2;

        double ce1;
        double se1;
        double ce2;
        double se2;

        double us_edge_k;
        double ds_edge_k;

        double us_edge_k_p;
        double ds_edge_k_p;

        double us_edge_k_x;
        double us_edge_k_y;

        double ds_edge_k_x;
        double ds_edge_k_y;

        Kokkos::complex<double> phase;
        Kokkos::complex<double> term;
    };

    struct PropSbend
    {
        Particles p;
        const SbendParams sp;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (sp.ledge)
            {
                // slot
                FF_algorithm::slot_unit( 
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                // edge:
                // 1. chef fixed angle (only kicks yp)
                // FF_algorithm::edge_unit(
                //        p(i,2), p(i,3), us_edge_k );
                //
                // 2. chef per-particle angle
                // FF_algorithm::edge_unit(
                //        p(i,2), p(i,1), p(i,3), dpop, us_edge_k_p);
                //
                // 3. ref particle angle (kicks both xp and yp)
                FF_algorithm::edge_unit(
                        p(i,2), p(i,1), p(i,3), sp.us_edge_k_x, sp.us_edge_k_y, 0); 
            }

            // bend
            FF_algorithm::bend_unit(
                    p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                    sp.dphi, sp.strength, sp.pref_b, sp.m_b, sp.ref_cdt, 
                    sp.phase, sp.term );

            if (sp.redge)
            {
                // edge
                // FF_algorithm::edge_unit(y, yp, ds_edge_k);
                // FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                FF_algorithm::edge_unit(
                        p(i,2), p(i,1), p(i,3), sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 

                // slot
                FF_algorithm::slot_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        sp.ce2, sp.se2, sp.pref_b, sp.m_b);
            }
        }
    };

    struct PropSbendCF
    {
        Particles p;
        const SbendParams sp;

        const double step_angle;
        const double step_ref_cdt;

        const Kokkos::complex<double> phase_e1;
        const Kokkos::complex<double> phase_e2;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (sp.ledge)
            {
                // slot
                FF_algorithm::slot_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                // edge
                //FF_algorithm::edge_unit(y, yp, us_edge_k);
                //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                FF_algorithm::edge_unit(
                        p(i,2), p(i,1), p(i,3), sp.us_edge_k_x, sp.us_edge_k_y, 0); 

                // bend edge (thin, but with face angle)
                FF_algorithm::bend_edge(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        sp.e1, phase_e1, sp.strength, sp.pref_b, sp.m_b);
            }

            // bend body
            FF_algorithm::bend_yoshida6< double, 
                                         FF_algorithm::thin_cf_kick_2<double>, 
                                         FF_algorithm::sbend_unit_phase,
                                         FF_algorithm::sbend_unit_term,
                                         FF_algorithm::sbend_dphi,
                                         2 > ( 
                  p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                  sp.pref_b, sp.m_b, step_ref_cdt,
                  step_angle, sp.step_kl,
                  sp.r0, sp.strength, sp.steps );

            if (sp.redge)
            {
                // bend edge (thin, but with face angle)
                FF_algorithm::bend_edge(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        sp.e2, phase_e2, sp.strength, sp.pref_b, sp.m_b);

                // edge
                //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                FF_algorithm::edge_unit(
                        p(i,2), p(i,1), p(i,3), sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 

                // slot
                FF_algorithm::slot_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        sp.ce2, sp.se2, sp.pref_b, sp.m_b);
            }
        }
    };

    void prop_reference(Reference_particle & ref_l, SbendParams & sp, double scale)
    {
        double pref_l = ref_l.get_momentum();
        double    m_l = ref_l.get_mass();

        // propagate the reference particle, and set the edge kick strength 
        // from the reference particle
        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        if (sp.ledge)
        {
            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce1, sp.se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.us_edge_k_x = sp.us_edge_k_p * (xp_l/zp_l);
            sp.us_edge_k_y = sp.us_edge_k_p * (yp_l/zp_l);

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.us_edge_k_x/scale, sp.us_edge_k_y/scale, 0);
        }

        FF_algorithm::bend_complete(
                x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
                sp.dphi, sp.strength, pref_l, m_l, 0.0/*ref cdt*/, 
                sp.phase, sp.term);

        if (sp.redge)
        {
            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.ds_edge_k_x = sp.ds_edge_k_p * (xp_l/zp_l);
            sp.ds_edge_k_y = sp.ds_edge_k_p * (yp_l/zp_l);

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            FF_algorithm::edge_unit(
                    y_l, xp_l, yp_l, sp.ds_edge_k_x/scale, sp.ds_edge_k_y/scale, 0);

            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce2, sp.se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        sp.ref_cdt = cdt_l;
    }

    void prop_reference_cf(
            Reference_particle & ref_l, SbendParams & sp, double scale,
            Kokkos::complex<double> const& phase_e1,
            Kokkos::complex<double> const& phase_e2,
            double const* step_kl )
    {
        double pref_l = ref_l.get_momentum();
        double    m_l = ref_l.get_mass();

        // propagate the reference particle, and set the edge kick strength 
        // from the reference particle
        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        if (sp.ledge)
        {
            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce1, sp.se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.us_edge_k_x = sp.us_edge_k_p * (xp_l/zp_l);
            sp.us_edge_k_y = sp.us_edge_k_p * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, us_edge_k/scale);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, us_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.us_edge_k_x/scale, sp.us_edge_k_y/scale, 0);

            // bend edge (thin)
            FF_algorithm::bend_edge(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.e1, phase_e1, sp.strength, pref_l, m_l);
        }

        FF_algorithm::bend_yoshida6< double, 
                                     FF_algorithm::thin_cf_kick_2<double>, 
                                     FF_algorithm::sbend_unit_phase,
                                     FF_algorithm::sbend_unit_term,
                                     FF_algorithm::sbend_dphi,
                                     2 >
            ( x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
              pref_l, m_l, 0.0 /* step ref_cdt */,
              sp.angle/sp.steps, step_kl,
              sp.r0, sp.strength, sp.steps);

        if (sp.redge)
        {
            // bend edge (thin)
            FF_algorithm::bend_edge(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.e2, phase_e2, sp.strength, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.ds_edge_k_x = sp.ds_edge_k_p * (xp_l/zp_l);
            sp.ds_edge_k_y = sp.ds_edge_k_p * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, ds_edge_k);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, ds_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.ds_edge_k_x/scale, sp.ds_edge_k_y/scale, 0);

            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce2, sp.se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        sp.ref_cdt = cdt_l;
    }
}

void FF_sbend::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF sbend on JetParticle not implemented");

#if 0
    double length = slice.get_right() - slice.get_left();
    double angle = slice.get_lattice_element().get_double_attribute("angle");

    double cos_angle = cos(angle);
    double sin_angle = sin(angle);

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

    sbend_unit(x, xp, y, yp, cdt, dpop,
               length, cos_angle, sin_angle,
               reference_momentum, m, reference_brho);
#endif
}

void FF_sbend::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    auto const& ele = slice.get_lattice_element();

    SbendParams sp;

    double      a = ele.get_double_attribute("angle");
    double      l = ele.get_double_attribute("l");

    sp.length = slice.get_right() - slice.get_left();
    sp.angle  = ( sp.length / l ) * a;
    sp.r0     = l / a;
    sp.steps  = steps;

    sp.ledge = slice.has_left_edge();
    sp.redge = slice.has_right_edge();

    sp.e1 = ele.get_double_attribute("e1", 0.0);
    sp.e2 = ele.get_double_attribute("e2", 0.0);

    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    int cf = 0;  // combined function
    double k_l[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // quad, sextupole, and octupole components
    k_l[0] = ele.get_double_attribute("k1", 0.0); k_l[1] = 0.0;
    k_l[2] = ele.get_double_attribute("k2", 0.0); k_l[3] = 0.0;
    k_l[4] = ele.get_double_attribute("k3", 0.0); k_l[5] = 0.0;

    if (k_l[0] != 0.0) cf = 1;
    if (k_l[2] != 0.0) cf = 2;
    if (k_l[4] != 0.0) cf = 3;

    double scale = ref_l.get_momentum() / 
        (ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]));

    double k_b[6] = { k_l[0] * scale, k_l[1] * scale,
                      k_l[2] * scale, k_l[3] * scale,
                      k_l[4] * scale, k_l[5] * scale };

    double usAngle = sp.e1;
    double dsAngle = -sp.e2;
    double usFaceAngle = sp.e1;
    double dsFaceAngle = sp.e2;

    if (!sp.redge)
    {
        dsAngle = 0;
        dsFaceAngle = 0;
    }

    if (!sp.ledge)
    {
        usAngle = 0;
        usFaceAngle = 0;
    }

    // lattice reference
    double pref_l = ref_l.get_momentum();
    double brho_l = pref_l / PH_CNV_brho_to_p;
    double    m_l = ref_l.get_mass();
    int  charge_l = ref_l.get_charge();

    // the particle dpop is with respect to this momentum which goes with the bunch
    double pref_b = ref_b.get_momentum();
    double brho_b = pref_b / PH_CNV_brho_to_p;
    double    m_b = bunch.get_mass();

    // common
    double psi = sp.angle - (usFaceAngle + dsFaceAngle);

    sp.strength = brho_l * a / l;
    sp.pref_b   = pref_b;
    sp.m_b      = m_b;

    sp.dphi  = -psi;
    sp.phase = std::exp( std::complex<double>(0.0, psi) );
    sp.term  = std::complex<double>(0.0, sp.length / sp.angle) *
               std::complex<double>(1.0 - cos(sp.angle), - sin(sp.angle)) *
               std::complex<double>(cos(dsFaceAngle), -sin(dsFaceAngle));

    sp.ce1 = cos(-sp.e1);
    sp.se1 = sin(-sp.e1);
    sp.ce2 = cos(-sp.e2);
    sp.se2 = sin(-sp.e2);

    // edge kick strength (scaled to particle)
    sp.us_edge_k =   ((charge_l > 0) ? 1.0 : -1.0) * sp.strength * tan(usAngle) / brho_b;
    sp.ds_edge_k = - ((charge_l > 0) ? 1.0 : -1.0) * sp.strength * tan(dsAngle) / brho_b;

    // edge kick (per) particle (full strength, angles are calculated at each particle)
    sp.us_edge_k_p =   ((charge_l > 0) ? 1.0 : -1.0) * sp.strength / brho_b;
    sp.ds_edge_k_p = - ((charge_l > 0) ? 1.0 : -1.0) * sp.strength / brho_b;

    // edge kick x/y
    sp.us_edge_k_x = 0.0;
    sp.us_edge_k_y = 0.0;

    sp.ds_edge_k_x = 0.0;
    sp.ds_edge_k_y = 0.0;


    if (cf == 0)
    {
        prop_reference(ref_l, sp, scale);

        int num = bunch.get_local_num(ParticleGroup::regular);
        auto parts = bunch.get_local_particles(ParticleGroup::regular);

        PropSbend sbend{parts, sp};
        Kokkos::parallel_for(num, sbend);
    }
    else
    {
        Kokkos::complex<double> phase_e1 = FF_algorithm::bend_edge_phase(sp.e1);
        Kokkos::complex<double> phase_e2 = FF_algorithm::bend_edge_phase(sp.e2);
        double step_length = sp.length / sp.steps;

        // step strength for reference
        for (int i=0; i<6; ++i) sp.step_kl[i] = k_l[i] * step_length;

        // propagate reference
        prop_reference_cf(ref_l, sp, scale, phase_e1, phase_e2, sp.step_kl);

        // step strength for bunch
        for (int i=0; i<6; ++i) sp.step_kl[i] = k_b[i] * step_length;

        // propagate bunch regular particles
        int num = bunch.get_local_num(ParticleGroup::regular);
        auto parts = bunch.get_local_particles(ParticleGroup::regular);

        PropSbendCF sbend{parts, sp, sp.angle/sp.steps, sp.ref_cdt/sp.steps, phase_e1, phase_e2};
        Kokkos::parallel_for(num, sbend);
    }

    bunch.get_reference_particle().increment_trajectory(sp.length);
}

