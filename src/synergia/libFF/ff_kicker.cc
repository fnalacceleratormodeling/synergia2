#include "ff_kicker.h"
#include "ff_algorithm.h"
#include "synergia/utils/simple_timer.h"

#include "ff_patterned_propagator.h"

void FF_kicker::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error(
            "libFF::hkicker for JetParticles has yet to be implemented");
}

void FF_kicker::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    scoped_simple_timer timer("libFF_kicker");

    auto const& elem = slice.get_lattice_element();
    const double length = slice.get_right() - slice.get_left();

    // hk and vk are the hk/vk under lattice reference momentum
    double    l = elem.get_double_attribute("l");
    double  hk0 = elem.get_double_attribute("hkick");
    double  vk0 = elem.get_double_attribute("vkick");
    double tilt = elem.get_double_attribute("tilt");

    // tilt
    double hk = cos(tilt)*hk0 - sin(tilt)*vk0;
    double vk = sin(tilt)*hk0 + cos(tilt)*vk0;

    //double k = elem.get_double_attribute(k_attr, 0.0);

    auto& ref_lattice = bunch.get_design_reference_particle();
    auto const& ref_bunch = bunch.get_reference_particle();

    double plattice = ref_lattice.get_momentum();
    double pbunch = ref_bunch.get_momentum();

    // scale is to scale the kick strength defined relative to the lattice momentum to
    // the scale of the bunch particles defined relative to the bunch momentum
    double scale = plattice/pbunch;

    // kick strength is defined as momentum change/reference momentum
    double b_hk = hk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;
    double b_vk = vk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;

    double k[2] = {hk, vk};
    double sk[2] = {b_hk, b_vk};

    using pp = FF_patterned_propagator<double, 1,
          FF_algorithm::thin_kicker_unit<double>>;

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference 
        // particle which is relative to the p-lattice.
        // also update the reference particle
        pp::get_reference_cdt_zero(ref_lattice, k);

        pp::apply_thin_kick(bunch, ParticleGroup::regular, sk);
        pp::apply_thin_kick(bunch, ParticleGroup::spectator, sk);
    }
    else
    {
        // simple drift-kick-drift scheme
        double simple_d = slice
            .get_lattice_element()
            .get_double_attribute("simple", 0.0);

        bool simple = fabs(simple_d) > 1e-16;

        // strength per unit length
        k[0] = k[0]/l;
        k[1] = k[1]/l;

        sk[0] = sk[0]/l;
        sk[1] = sk[1]/l;

        double pref = bunch.get_reference_particle().get_momentum();
        double mass = bunch.get_mass();

        if (simple)
        {
            // use un-scaled k_pul
            double ref_cdt = pp::get_reference_cdt_simple(
                    ref_lattice, length, k);

            pp::apply_simple_kick(bunch, ParticleGroup::regular, 
                    pref, mass, ref_cdt, length, sk);

            pp::apply_simple_kick(bunch, ParticleGroup::spectator, 
                    pref, mass, ref_cdt, length, sk);
        }
        else
        {
            // yoshida steps
            steps = (int)slice
                .get_lattice_element()
                .get_double_attribute("yoshida_steps", 4.0);

            // use un-scaled k_pul
            double ref_cdt = pp::get_reference_cdt_yoshida(
                    ref_lattice, length, k, steps, true);

            pp::apply_yoshida_kick(bunch, ParticleGroup::regular,
                    pref, mass, ref_cdt, length, sk, steps);

            pp::apply_yoshida_kick(bunch, ParticleGroup::spectator,
                    pref, mass, ref_cdt, length, sk, steps);
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
}

