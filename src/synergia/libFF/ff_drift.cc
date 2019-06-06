#include "ff_drift.h"

#include "synergia/libFF/ff_algorithm.h"
//#include "synergia/lattice/chef_utils.h"
//#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

namespace
{
    struct PropDrift
    {
        Particles p;
        double l, ref_p, m, ref_t;

        PropDrift(Particles p_, double l, double ref_p, double m, double ref_t)
            : p(p_), l(l), ref_p(ref_p), m(m), ref_t(ref_t) { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            FF_algorithm::drift_unit(
                    p(i, 0), p(i, 1), p(i, 2),
                    p(i, 3), p(i, 4), p(i, 5),
                    l, ref_p, m, 0.0);
        }
    };

    double get_reference_cdt(double length, Reference_particle & ref)
    {
        double x(ref.get_state()[Bunch::x]);
        double xp(ref.get_state()[Bunch::xp]);
        double y(ref.get_state()[Bunch::y]);
        double yp(ref.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(ref.get_state()[Bunch::dpop]);
        double ref_p = ref.get_momentum();
        double m = ref.get_mass();

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, m, 0.0);

        // propagate and update the bunch design reference particle state
        ref.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}

void FF_drift::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a drift element is yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

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
    double reference_cdt = get_reference_cdt(length, reference_particle);

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
#endif
}

void FF_drift::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const double    mass = bunch.get_mass();

    const int local_num   = bunch.get_local_num();
    const int local_s_num = bunch.get_local_spectator_num();

    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();
    const double ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    const double ref_cdt = get_reference_cdt(length, ref_l);

    auto parts = bunch.get_local_particles();
    auto hparts = bunch.get_host_particles();

    Kokkos::parallel_for(
            local_num, PropDrift(parts, length, ref_p, mass, ref_cdt) );

    bunch.checkout_particles();

    for (int p=0; p<local_num; ++p)
    {
        std::cout 
            << hparts(p, 0) << ", "
            << hparts(p, 1) << ", "
            << hparts(p, 2) << ", "
            << hparts(p, 3) << ", "
            << hparts(p, 4) << ", "
            << hparts(p, 5) << "\n";
    }

    bunch.get_reference_particle().increment_trajectory(length);

#if 0
    // real particles
    {
        bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

        const int gsvsize = GSVector::size();
        const int num_blocks = local_num / gsvsize;
        const int block_last = num_blocks * gsvsize;

        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize)
        {
            GSVector x(&xa[part]);
            GSVector xp(&xpa[part]);
            GSVector y(&ya[part]);
            GSVector yp(&ypa[part]);
            GSVector cdt(&cdta[part]);
            GSVector dpop(&dpopa[part]);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

            x.store(&xa[part]);
            y.store(&ya[part]);
            cdt.store(&cdta[part]);
        }

        for (int part = block_last; part < local_num; ++part)
        {
            double x(xa[part]);
            double xp(xpa[part]);
            double y(ya[part]);
            double yp(ypa[part]);
            double cdt(cdta[part]);
            double dpop(dpopa[part]);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

            xa[part] = x;
            ya[part] = y;
            cdta[part] = cdt;
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }

    // spectators
    {
        bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

        const int gsvsize = GSVector::size();
        const int num_blocks = local_s_num / gsvsize;
        const int block_last = num_blocks * gsvsize;

        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize)
        {
            GSVector x(&xa[part]);
            GSVector xp(&xpa[part]);
            GSVector y(&ya[part]);
            GSVector yp(&ypa[part]);
            GSVector cdt(&cdta[part]);
            GSVector dpop(&dpopa[part]);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

            x.store(&xa[part]);
            y.store(&ya[part]);
            cdt.store(&cdta[part]);
        }

        for (int part = block_last; part < local_s_num; ++part)
        {
            double x(xa[part]);
            double xp(xpa[part]);
            double y(ya[part]);
            double yp(ypa[part]);
            double cdt(cdta[part]);
            double dpop(dpopa[part]);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

            xa[part] = x;
            ya[part] = y;
            cdta[part] = cdt;
        }
    }
#endif
}

