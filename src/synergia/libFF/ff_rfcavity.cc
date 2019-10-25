#include "ff_rfcavity.h"
#include "ff_algorithm.h"

#include "synergia/foundation/math_constants.h"
#include "synergia/utils/simple_timer.h"

namespace
{
    struct RFCavityParams
    {
        int nh;
        double mhp[12];

        double length;
        double str;
        double phi_s;
        double w_rf;

        double pref_b;
        double m_b;
        double new_pref_b;

        double ref_cdt_1;
        double ref_cdt_2;
    };

    struct PropRFCavity
    {
        Particles p;
        ConstParticleMasks masks;
        const RFCavityParams rp;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                FF_algorithm::drift_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        0.5 * rp.length, rp.pref_b, rp.m_b, rp.ref_cdt_1);

                FF_algorithm::thin_rfcavity_unit(
                        p(i,1), p(i,3), p(i,4), p(i,5),
                        rp.w_rf, rp.str, rp.phi_s, rp.m_b, rp.pref_b, rp.new_pref_b, 
                        rp.mhp, rp.nh);

                FF_algorithm::drift_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        0.5 * rp.length, rp.pref_b, rp.m_b, rp.ref_cdt_2);
            }
        }
    };


}

#if 0
double get_reference_cdt(double length, Reference_particle & reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);
    double reference_momentum = reference_particle.get_momentum();
    double m = reference_particle.get_mass();

    // By definition, the reference particle receives no kick in the RF cavity because it is
    // perfectly in sync with the RF waveform.
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m, 0.0);
    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}
#endif

void FF_rfcavity::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    scoped_simple_timer timer("libFF_rfcavity");

    RFCavityParams rp;

    rp.length = slice.get_right() - slice.get_left();
    Lattice_element const & elm = slice.get_lattice_element();

    int    harmonic_number = elm.get_double_attribute("harmon", -1.0);
    double            volt = elm.get_double_attribute("volt", 0.0);
    double             lag = elm.get_double_attribute("lag", 0.0);
    double           shunt = elm.get_double_attribute("shunt", 0.0);
    // freq is the synchronous frequency of the cavity in MHz
    double            freq = elm.get_double_attribute("freq", -1.0);
    // delta_freq is the frequency offset from synchronous in MHz like freq
    double      delta_freq = elm.get_double_attribute("delta_freq", 0.0);

    // harmonics, 
    // mhp[h*3]   = harmonic multiple
    // mhp[h*3+1] = relative strength
    // mhp[h*3+2] = phase shift
    // nh = number of harmonics
    rp.nh = 1;
    rp.mhp[0] = 1.0; 
    rp.mhp[1] = 1.0; 
    rp.mhp[2] = 0.0;

    // harmonic multiple = 2
    if (elm.has_double_attribute("h2_str"))
    {
        rp.mhp[rp.nh*3+0] = 2; 
        rp.mhp[rp.nh*3+1] = elm.get_double_attribute("h2_str", 0.0);
        rp.mhp[rp.nh*3+2] = elm.get_double_attribute("h2_phase", 0.0);
        ++rp.nh;
    }

    // harmonic multiple = 3
    if (elm.has_double_attribute("h3_str"))
    {
        rp.mhp[rp.nh*3+0] = 3; 
        rp.mhp[rp.nh*3+1] = elm.get_double_attribute("h3_str", 0.0);
        rp.mhp[rp.nh*3+2] = elm.get_double_attribute("h3_phase", 0.0);
        ++rp.nh;
    }

    // harmonic multiple = 4
    if (elm.has_double_attribute("h4_str"))
    {
        rp.mhp[rp.nh*3+0] = 4; 
        rp.mhp[rp.nh*3+1] = elm.get_double_attribute("h4_str", 0.0);
        rp.mhp[rp.nh*3+2] = elm.get_double_attribute("h4_phase", 0.0);
        ++rp.nh;
    }

    rp.str = volt * 1.0e-3;

    // keep lag within the range of [0, 1).
    while (lag < 0.0)  { lag += 1.0; }
    while (lag >= 1.0) { lag -= 1.0; }

    //elm.set_double_attribute("lag", lag);
    rp.phi_s = 2.0 * mconstants::pi * lag;
    rp.w_rf  = 2.0 * mconstants::pi * (freq + delta_freq) * 1.0e6;

    Reference_particle & ref_l = bunch.get_design_reference_particle();
    Reference_particle & ref_b = bunch.get_reference_particle();

    // The bunch particles momentum is with respect to the bunch reference particle
    rp.pref_b = ref_b.get_momentum();
    rp.m_b = bunch.get_mass();

    rp.new_pref_b = FF_algorithm::thin_rfcavity_pnew(
            rp.pref_b, rp.m_b, rp.str, rp.phi_s);

    // reference_cdt uses the lattice reference particle
    // double reference_cdt = get_reference_cdt(length, ref_l);
    double ref_l_x    = ref_l.get_state()[Bunch::x];
    double ref_l_xp   = ref_l.get_state()[Bunch::xp];
    double ref_l_y    = ref_l.get_state()[Bunch::y];
    double ref_l_yp   = ref_l.get_state()[Bunch::yp];
    double ref_l_cdt  = 0.0;
    double ref_l_dpop = ref_l.get_state()[Bunch::dpop];

    double ref_l_p = ref_l.get_momentum();
    double ref_l_m = ref_l.get_mass();

    // double new_ref_l_p = 
    //     FF_algorithm::thin_rfcavity_pnew(ref_l_p, ref_l_m, str, phi_s);

    // first half drift
    FF_algorithm::drift_unit(
            ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop,
            0.5 * rp.length, ref_l_p, ref_l_m, 0.0);

    double total_ref_cdt = ref_l_cdt;
    rp.ref_cdt_1 = ref_l_cdt;
    ref_l_cdt = 0.0;

    // do not give it the new_ref_l_p because the xp and yp dont get scaled, the momentum 
    // of the lattice reference particle remains unchanged, only the dpop of the state
    // has been changed
    FF_algorithm::thin_rfcavity_unit(
            ref_l_xp, ref_l_yp, ref_l_cdt, ref_l_dpop,
            rp.w_rf, rp.str, rp.phi_s, ref_l_m, ref_l_p, ref_l_p, 
            rp.mhp, rp.nh );

    // second half drift
    FF_algorithm::drift_unit(
            ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop,
            0.5 * rp.length, ref_l_p, ref_l_m, 0.0);

    total_ref_cdt += ref_l_cdt;
    rp.ref_cdt_2 = ref_l_cdt;

    // save the state
    ref_l.set_state(ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, total_ref_cdt, ref_l_dpop);

    // bunch particles
    int num = bunch.get_local_num_slots(ParticleGroup::regular);
    auto parts = bunch.get_local_particles(ParticleGroup::regular);
    auto masks = bunch.get_local_particles_masks(ParticleGroup::regular);

    PropRFCavity rfcavity{parts, masks, rp};
    Kokkos::parallel_for(num, rfcavity);

    // bunch spectator particles
    // TODO: ...

    // updated four momentum
    Four_momentum fm = ref_b.get_four_momentum();
    fm.set_momentum(rp.new_pref_b);

    // update the bunch reference particle with the updated ref_p
    ref_b.set_four_momentum(fm);
    ref_b.increment_trajectory(rp.length);
}

