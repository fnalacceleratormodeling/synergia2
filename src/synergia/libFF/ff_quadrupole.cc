#include "ff_quadrupole.h"

#include "synergia/libFF/ff_algorithm.h"
//#include "synergia/lattice/chef_utils.h"
//#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

namespace
{
    struct PropQuadThin
    {
        Particles p;
        double k[2];
        double xoff, yoff;

        PropQuadThin(Particles p_, double k0, double k1, double xoff, double yoff)
            : p(p_), k{k0, k1}, xoff(xoff), yoff(yoff) { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            double x = p(i, 0) - xoff;
            double y = p(i, 2) - yoff;

            FF_algorithm::thin_quadrupole_unit(
                    x, p(i, 1), y, p(i, 3), k);
        }
    };

    struct PropQuad
    {
        Particles p;
        int steps;
        double xoff, yoff;
        double ref_p, ref_m, step_ref_t, step_l, step_k[2];

        PropQuad( Particles p_, 
                  int steps,
                  double xoff,
                  double yoff,
                  double ref_p, 
                  double ref_m, 
                  double ref_t, 
                  double length,
                  double k0,
                  double k1 )
            : p(p_)
            , steps(steps)
            , xoff(xoff)
            , yoff(yoff)
            , ref_p(ref_p)
            , ref_m(ref_m)
            , step_ref_t(ref_t/steps)
            , step_l(length/steps) 
            , step_k{k0*step_l, k1*step_l}
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            p(i, 0) -= xoff;
            p(i, 2) -= yoff;

            FF_algorithm::yoshida6<double, FF_algorithm::thin_quadrupole_unit<double>, 1>
                (p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                 ref_p, ref_m, step_ref_t, step_l, step_k, steps);

            p(i, 0) += xoff;
            p(i, 2) += yoff;
        }
    };


    double get_reference_cdt(double length, int steps, double * k, Reference_particle &ref)
    {
        if (length == 0)
        {
            ref.set_state_cdt(0.0);
            return 0.0;
        }
        else
        {
            double ref_p = ref.get_momentum();
            double ref_m = ref.get_mass();

            // steps comes from base class, set in apply
            double step_length = length/steps;
            double step_strength[2] = { k[0]*step_length, k[1]*step_length };

            double x(ref.get_state()[Bunch::x]);
            double xp(ref.get_state()[Bunch::xp]);
            double y(ref.get_state()[Bunch::y]);
            double yp(ref.get_state()[Bunch::yp]);
            double cdt(0.0);
            double dpop(ref.get_state()[Bunch::dpop]);

            FF_algorithm::yoshida4<double, FF_algorithm::thin_quadrupole_unit<double>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      ref_p, ref_m, 0.0,
                      step_length, step_strength, steps );

            // propagate and update the lattice reference particle state
            ref.set_state(x, xp, y, yp, cdt, dpop);

            return cdt;
        }
    }
}

void FF_quadrupole::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a quadrupole element is yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k[1] = slice.get_lattice_element().get_double_attribute("k1s", 0.0);

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
    double m = jet_particle.Mass();

    Particle chef_particle(jet_particle);
    Reference_particle reference_particle(
                chef_particle_to_reference_particle(chef_particle));
    double reference_cdt = get_reference_cdt(length, k, reference_particle);
    double step_reference_cdt = reference_cdt/steps;
    double step_length = length/steps;
    double step_strength[2] = { k[0]*step_length, k[1]*step_length };
    double kl[2] = { k[0]*length, k[1]*length };

    if (length == 0.0) {
        FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kl);
    } else {
        FF_algorithm::yoshida4<TJet<double>, FF_algorithm::thin_quadrupole_unit<TJet<double> >, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  step_reference_cdt,
                  step_length, step_strength, steps );
    }
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
#endif
}

void FF_quadrupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    auto const& ele = slice.get_lattice_element();

    // length
    double length = slice.get_right() - slice.get_left();

    // strength
    double k[2];
    k[0] = ele.get_double_attribute("k1", 0.0);
    k[1] = ele.get_double_attribute("k1s", 0.0);

    // offsets
    const double xoff = ele.get_double_attribute("hoffset", 0.0);
    const double yoff = ele.get_double_attribute("voffset", 0.0);

    // tilting
    double tilt = ele.get_double_attribute("tilt", 0.0);
    if (tilt != 0.0)
    {
        std::complex<double> ck2(k[0], -k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -2.0*tilt));
        k[0] = ck2.real();
        k[1] = ck2.imag();
    }

    // scaling
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
    double brho_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

    double scale = brho_l / brho_b;

    k[0] *= scale;
    k[1] *= scale;

    if (close_to_zero(length))
    {
        // TODO: should move this to the get_reference_cdt() method
        // propagate the bunch design reference particle
        double x  = ref_l.get_state()[Bunch::x];
        double xp = ref_l.get_state()[Bunch::xp];
        double y  = ref_l.get_state()[Bunch::y];
        double yp = ref_l.get_state()[Bunch::yp];
        double dpop = ref_l.get_state()[Bunch::dpop];

        x -= xoff;
        y -= yoff;

        FF_algorithm::thin_quadrupole_unit(x, xp,  y, yp, k);

        x += xoff;
        y += yoff;

        ref_l.set_state(x, xp, y, yp, 0.0, dpop);

        // propagate the bunch particles
        int num = bunch.get_local_num(ParticleGroup::regular);
        auto parts = bunch.get_local_particles(ParticleGroup::regular);
        Kokkos::parallel_for(num, PropQuadThin(parts, k[0], k[1], xoff, yoff) );

        // TODO: spectator particles
        // ...
    }
    else
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // params
        double ref_p = ref_b.get_momentum();
        double ref_m = ref_b.get_mass();
        double ref_t = get_reference_cdt(length, steps, k, ref_l);

        // bunch particles
        int num = bunch.get_local_num(ParticleGroup::regular);
        auto parts = bunch.get_local_particles(ParticleGroup::regular);
        PropQuad pq( parts, steps, xoff, yoff, ref_p, ref_m, ref_t, length, k[0], k[1] );
        Kokkos::parallel_for(num, pq);

        // TODO: spectator particles
        // ...

        // advance the ref_part
        bunch.get_reference_particle().increment_trajectory(length);
    }
}

