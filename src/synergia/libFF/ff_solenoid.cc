#include "ff_solenoid.h"
#include "ff_drift.h"
#include "synergia/libFF/ff_algorithm.h"


namespace
{
    template<void(KF)(double const&, double&, double const&, double&, double)>
    struct solenoid_edge_kicker
    {
        Particles p;
        ParticleMasks m;

        double kse;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { if (m(i)) KF(p(i,0), p(i,1), p(i,2), p(i,3), kse); }
    };

    template<void(KF)(double const&, double&, double const&, double&, double)>
    void apply_edge_kick(Bunch& bunch, ParticleGroup pg, double kse)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        solenoid_edge_kicker<KF> sk{parts, masks, kse};
        Kokkos::parallel_for(bunch.size(pg), sk);
    }

    struct solenoid_unit_kicker
    {
        Particles p;
        ParticleMasks m;

        double ksl;
        double ks;
        double length;
        double ref_p;
        double mass;
        double ref_cdt;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if (m(i)) 
                FF_algorithm::solenoid_unit(
                        p(i,0), p(i,1), p(i,2), 
                        p(i,3), p(i,4), p(i,5),
                        ksl, ks, length, ref_p, mass, ref_cdt); 
        }
    };

    void apply_solenoid_unit(Bunch& bunch, ParticleGroup pg,
            double ksl, double ks, double length, double ref_p,
            double mass, double ref_cdt)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        solenoid_unit_kicker sk{
            parts, masks, ksl, ks, length, 
            ref_p, mass, ref_cdt
        };

        Kokkos::parallel_for(bunch.size(pg), sk);
    }


    double get_reference_cdt_solenoid(
            double length, Reference_particle & reference_particle,
            bool in_edge, bool out_edge, double ks, double kse, double ksl)
    {
        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double ref_p = reference_particle.get_momentum()*(1+dpop);
        double mass  = reference_particle.get_mass();

        // in edge
        if (in_edge)
        {
            FF_algorithm::solenoid_in_edge_kick(x, xp, y, yp, kse);
        }

        // body
        FF_algorithm::solenoid_unit( x, xp, y, yp, cdt, dpop,
                ksl, ks, length, ref_p, mass, 0.0);

        // out edge
        if (out_edge)
        {
            FF_algorithm::solenoid_out_edge_kick(x, xp, y, yp, kse);
        }

        // propagate and update the bunch design reference particle state
        reference_particle.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}



void FF_solenoid::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error(
            "Propagate JetParticle through a solenoid element is yet to be implemented");
}

void FF_solenoid::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const double    mass = bunch.get_mass();

    // in edge, out edge
    const bool has_in_edge  = slice.has_left_edge();
    const bool has_out_edge = slice.has_right_edge();

    // ks here is the field
    double ks = slice.get_lattice_element().get_double_attribute("ks");

    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();
    const double ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    // if ks=0 it is effectively a drift
    if (fabs(ks) < 1e-12)
    {
        FF_drift drift;
        drift.apply(slice, bunch);

        return;
    }

    // scaling
    double brho_l = ref_l.get_momentum() 
        * (1.0 + ref_l.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

    double brho_b = ref_b.get_momentum() 
        * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

    double scale = brho_l / brho_b;

            ks = ks * scale;
    double ksl = ks * length;
    double kse = ks * 0.5;

    // reference cdt
    const double ref_cdt = get_reference_cdt_solenoid(length, ref_l, 
            has_in_edge, has_out_edge, ks, kse, ksl);

    // in-edge
    if (has_in_edge)
    {
        apply_edge_kick<FF_algorithm::solenoid_in_edge_kick>(
                bunch, ParticleGroup::regular, kse);

        apply_edge_kick<FF_algorithm::solenoid_in_edge_kick>(
                bunch, ParticleGroup::spectator, kse);
    }

    // body
    apply_solenoid_unit(bunch, ParticleGroup::regular,
            ksl, ks, length, ref_p, mass, ref_cdt);

    apply_solenoid_unit(bunch, ParticleGroup::spectator,
            ksl, ks, length, ref_p, mass, ref_cdt);

    // out-edge
    if (has_out_edge)
    {
        apply_edge_kick<FF_algorithm::solenoid_out_edge_kick>(
                bunch, ParticleGroup::regular, kse);

        apply_edge_kick<FF_algorithm::solenoid_out_edge_kick>(
                bunch, ParticleGroup::spectator, kse);
    }

    // trajectory
    bunch.get_reference_particle().increment_trajectory(length);
}

