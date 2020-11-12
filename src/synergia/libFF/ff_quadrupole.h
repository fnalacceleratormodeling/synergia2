#ifndef FF_QUADRUPOLE_H
#define FF_QUADRUPOLE_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"

namespace quad_impl
{
    template<class BP>
    struct PropQuadThin
    {
        typename BP::parts_t p;
        ConstParticleMasks masks;
        double k[2];
        double xoff, yoff;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                auto x = p(i, 0) - xoff;
                auto y = p(i, 2) - yoff;

                FF_algorithm::thin_quadrupole_unit(
                        x, p(i, 1), y, p(i, 3), k);
            }
        }
    };

    template<class BP>
    struct PropQuad
    {
        typename BP::parts_t p;
        ConstParticleMasks masks;
        int steps;
        double xoff, yoff;
        double ref_p, ref_m, step_ref_t, step_l, step_k[2];

        PropQuad( BP & bp,
                  int steps,
                  double xoff,
                  double yoff,
                  double ref_p, 
                  double ref_m, 
                  double ref_t, 
                  double length,
                  double k0,
                  double k1 )
            : p(bp.parts)
            , masks(bp.masks)
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
            if (masks(i))
            {
                p(i, 0) -= xoff;
                p(i, 2) -= yoff;

                FF_algorithm::yoshida6<typename BP::part_t, 
                    FF_algorithm::thin_quadrupole_unit<typename BP::part_t>, 1> (
                            p(i,0), p(i,1), p(i,2), 
                            p(i,3), p(i,4), p(i,5), 
                            ref_p, ref_m, step_ref_t, 
                            step_l, step_k, steps );

                p(i, 0) += xoff;
                p(i, 2) += yoff;
            }
        }
    };

    template<class BP>
    struct PropQuadSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        ConstParticleMasks masks;
        int steps;
        gsv_t xoff, yoff;
        double ref_p, ref_m, step_ref_t, step_l, step_k[2];

        PropQuadSimd( BP & bp,
                  int steps,
                  double xoff,
                  double yoff,
                  double ref_p, 
                  double ref_m, 
                  double ref_t, 
                  double length,
                  double k0,
                  double k1 )
            : p(bp.parts)
            , masks(bp.masks)
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
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                p0 = p0 - xoff;
                p2 = p2 - yoff;

                FF_algorithm::yoshida6<gsv_t, 
                    FF_algorithm::thin_quadrupole_unit<gsv_t>, 1>(
                            p0, p1, p2, p3, p4, p5,
                            ref_p, ref_m, step_ref_t, 
                            step_l, step_k, steps );

                p0 = p0 + xoff;
                p2 = p2 + yoff;

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
            }
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

            FF_algorithm::yoshida6<double, 
                FF_algorithm::thin_quadrupole_unit<double>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      ref_p, ref_m, 0.0,
                      step_length, step_strength, steps );

            // propagate and update the lattice reference particle state
            ref.set_state(x, xp, y, yp, cdt, dpop);

            return cdt;
        }
    }
}


namespace FF_quadrupole
{
    template<class BunchT>
    inline void apply(Lattice_element_slice const& slice, BunchT & bunch)
    {
        using namespace quad_impl;

        scoped_simple_timer timer("libFF_quad");

        // element
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
            std::complex<double> ck2(k[0], k[1]);
            ck2 = ck2 * exp(std::complex<double>(0.0, -2.0*tilt));
            k[0] = ck2.real();
            k[1] = ck2.imag();
        }

        // scaling
        Reference_particle       & ref_l = bunch.get_design_reference_particle();
        Reference_particle const & ref_b = bunch.get_reference_particle();

        double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
        double brho_b = ref_b.get_momentum()
                        * (1.0 + ref_b.get_state()[Bunch::dpop])
                        / ref_l.get_charge();  // GV/c

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
            auto bp = bunch.get_bunch_particles(ParticleGroup::regular);

            PropQuadThin<typename BunchT::bp_t> pqt{ bp.parts, bp.masks, 
                {k[0], k[1]}, xoff, yoff };
            Kokkos::parallel_for(bp.size(), pqt);

            // TODO: spectator particles
            // ...
        }
        else
        {
            // yoshida steps
            int steps = (int)slice.get_lattice_element()
                                  .get_double_attribute("yoshida_steps", 4.0);

            // params
            double ref_p = ref_b.get_momentum();
            double ref_m = ref_b.get_mass();
            double ref_t = get_reference_cdt(length, steps, k, ref_l);

            // bunch particles
            auto apply = [&](ParticleGroup pg) {
                auto bp = bunch.get_bunch_particles(pg);
                if (!bp.size()) return;

#if LIBFF_USE_GSV
                PropQuadSimd<typename BunchT::bp_t> pq(bp, steps, xoff, yoff, 
                        ref_p, ref_m, ref_t, length, k[0], k[1]);
                Kokkos::parallel_for(bp.size_in_gsv(), pq);
#else
                PropQuad<typename BunchT::bp_t> pq(bp, steps, xoff, yoff, 
                        ref_p, ref_m, ref_t, length, k[0], k[1]);
                Kokkos::parallel_for(bp.size(), pq);
#endif
            };

            // apply
            apply(ParticleGroup::regular);
            apply(ParticleGroup::spectator);

            // advance the ref_part
            bunch.get_reference_particle().increment_trajectory(length);
        }

        Kokkos::fence();
    }
}

#endif // FF_QUADRUPOLE_H
