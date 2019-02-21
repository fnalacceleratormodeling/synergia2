#include "ff_rfcavity.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/foundation/math_constants.h"

FF_rfcavity::FF_rfcavity()
{

}

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


void FF_rfcavity::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("propagate jet particle for rf cavity has yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();
    double freq = 0.0;

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

    rfcavity_unit(x, xp, y, yp, cdt, dpop,
                  length, freq,
                  reference_momentum, m, reference_brho);
#endif
}

void FF_rfcavity::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();
    Lattice_element const & elm = slice.get_lattice_element();

    int    harmonic_number = elm.get_double_attribute("harmon", -1.0);
    double            volt = elm.get_double_attribute("volt", 0.0);
    double             lag = elm.get_double_attribute("lag", 0.0);
    double           shunt = elm.get_double_attribute("shunt", 0.0);
    // freq is the synchronous frequency of the cavity in MHz
    double            freq = elm.get_double_attribute("freq", -1.0);
    // delta_freq is the frequency offset from synchronous in MHz like freq
    double            delta_freq = elm.get_double_attribute("delta_freq", 0.0);

    // harmonics, 
    // mhp[h*3]   = harmonic multiple
    // mhp[h*3+1] = relative strength
    // mhp[h*3+2] = phase shift
    // nh = number of harmonics
    double mhp[30]; int nh = 1;
    mhp[0] = 1.0; mhp[1] = 1.0; mhp[2] = 0.0;

    // harmonic multiple = 2
    if (elm.has_double_attribute("h2_str"))
    {
        mhp[nh*3+0] = 2; 
        mhp[nh*3+1] = elm.get_double_attribute("h2_str", 0.0);
        mhp[nh*3+2] = elm.get_double_attribute("h2_phase", 0.0);
        ++nh;
    }

    // harmonic multiple = 3
    if (elm.has_double_attribute("h3_str"))
    {
        mhp[nh*3+0] = 3; 
        mhp[nh*3+1] = elm.get_double_attribute("h3_str", 0.0);
        mhp[nh*3+2] = elm.get_double_attribute("h3_phase", 0.0);
        ++nh;
    }

    // harmonic multiple = 4
    if (elm.has_double_attribute("h4_str"))
    {
        mhp[nh*3+0] = 4; 
        mhp[nh*3+1] = elm.get_double_attribute("h4_str", 0.0);
        mhp[nh*3+2] = elm.get_double_attribute("h4_phase", 0.0);
        ++nh;
    }

    double   str = volt * 1.0e-3;

    // keep lag within the range of [0, 1).
    while (lag < 0.0)  { lag += 1.0; }
    while (lag >= 1.0) { lag -= 1.0; }

    //elm.set_double_attribute("lag", lag);
    double phi_s = 2.0 * mconstants::pi * lag;
    double  w_rf = 2.0 * mconstants::pi * (freq + delta_freq) * 1.0e6;

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    Reference_particle & ref_l = bunch.get_design_reference_particle();
    Reference_particle & ref_b = bunch.get_reference_particle();

    // The bunch particles momentum is with respect to the bunch reference particle
    double reference_momentum = ref_b.get_momentum();
    double m = bunch.get_mass();
    double new_ref_p = FF_algorithm::thin_rfcavity_pnew(reference_momentum, m, str, phi_s);

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

    double new_ref_l_p = FF_algorithm::thin_rfcavity_pnew(ref_l_p, ref_l_m, str, phi_s);

    // first half drift
    FF_algorithm::drift_unit(ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop,
            0.5 * length, ref_l_p, ref_l_m, 0.0);

    // do not give it the new_ref_l_p because the xp and yp dont get scaled, the momentum 
    // of the lattice reference particle remains unchanged, only the dpop of the state
    // has been changed
    FF_algorithm::thin_rfcavity_unit(ref_l_xp, ref_l_yp, ref_l_cdt, ref_l_dpop,
            w_rf, str, phi_s, ref_l_m, ref_l_p, ref_l_p, mhp, nh);

    // second half drift
    FF_algorithm::drift_unit(ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop,
            0.5 * length, ref_l_p, ref_l_m, 0.0);

    // save the state
    ref_l.set_state(ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop);

    // and finally the reference time
    double reference_cdt = ref_l_cdt;

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

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop,
                    0.5 * length, reference_momentum, m, 0.5 * reference_cdt);

            FF_algorithm::thin_rfcavity_unit(xp, yp, cdt, dpop,
                    w_rf, str, phi_s, m, reference_momentum, new_ref_p, mhp, nh);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop,
                    0.5 * length, reference_momentum, m, 0.5 * reference_cdt);

            particles[part][Bunch::x]    = x;
            particles[part][Bunch::xp]   = xp;
            particles[part][Bunch::y]    = y;
            particles[part][Bunch::yp]   = yp;
            particles[part][Bunch::cdt]  = cdt;
            particles[part][Bunch::dpop] = dpop;
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

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop,
                    0.5 * length, reference_momentum, m, 0.5 * reference_cdt);

            FF_algorithm::thin_rfcavity_unit(xp, yp, cdt, dpop,
                    w_rf, str, phi_s, m, reference_momentum, new_ref_p, mhp, nh);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop,
                    0.5 * length, reference_momentum, m, 0.5 * reference_cdt);

            s_particles[part][Bunch::x]    = x;
            s_particles[part][Bunch::xp]   = xp;
            s_particles[part][Bunch::y]    = y;
            s_particles[part][Bunch::yp]   = yp;
            s_particles[part][Bunch::cdt]  = cdt;
            s_particles[part][Bunch::dpop] = dpop;
        }
    }

    // updated four momentum
    Four_momentum fm = ref_b.get_four_momentum();
    fm.set_momentum(new_ref_p);

    // update the bunch reference particle with the updated ref_p
    ref_b.set_four_momentum(fm);
    ref_b.increment_trajectory(length);
}

template<class Archive>
    void
    FF_rfcavity::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_rfcavity::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_rfcavity::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_rfcavity::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_rfcavity::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_rfcavity::~FF_rfcavity()
{

}

