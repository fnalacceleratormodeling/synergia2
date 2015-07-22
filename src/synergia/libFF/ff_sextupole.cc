#include "ff_sextupole.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"

const int FF_sextupole::drifts_per_step = 4; // determined by algorithm in thick_sextupole unit

FF_sextupole::FF_sextupole()
{

}

double FF_sextupole::get_reference_cdt(double length, double * k,
                                       Reference_particle &reference_particle) {
    double reference_cdt;
    if (length == 0) {
        reference_cdt = 0.0;
    } else {
        double reference_momentum = reference_particle.get_momentum();
        double m = reference_particle.get_mass();
        double step_length = length/steps;
        double step_strength[2] = { k[0]*step_length, k[1]*step_length };

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(reference_particle.get_state()[Bunch::cdt]);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double cdt_orig = cdt;
        FF_algorithm::yoshida4<double, FF_sextupole::thin_sextupole_unit<double>, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  0.0,
                  step_length, step_strength, steps );
        reference_cdt = cdt - cdt_orig;
    }
    return reference_cdt;
}

void FF_sextupole::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    double length = slice.get_right() - slice.get_left();

    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k2");
    k[1] = slice.get_lattice_element().get_double_attribute("k2s");

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
    double substep_reference_cdt = reference_cdt/steps/drifts_per_step;
    double step_length = length/steps;
    double step_strength[2] = { k[0]*step_length, k[1]*step_length };
    double kl[2] = { k[0]*length, k[1]*length };

    if (length == 0.0) {
        thin_sextupole_unit(x, xp, y, yp, kl);
    } else {
        FF_algorithm::yoshida4<TJet<double>, FF_sextupole::thin_sextupole_unit<TJet<double> >, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  substep_reference_cdt,
                  step_length, step_strength, steps );
    }
    FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
}

void FF_sextupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
     double length = slice.get_right() - slice.get_left();

     double k[2];
     k[0] = slice.get_lattice_element().get_double_attribute("k2");
     k[1] = slice.get_lattice_element().get_double_attribute("k2s");

     int local_num = bunch.get_local_num();
     MArray2d_ref particles = bunch.get_local_particles();

     if (length == 0.0) {
         for (int part = 0; part < local_num; ++part) {
             double x(particles[part][Bunch::x]);
             double xp(particles[part][Bunch::xp]);
             double y(particles[part][Bunch::y]);
             double yp(particles[part][Bunch::yp]);

             thin_sextupole_unit(x, xp, y, yp, k);

             particles[part][Bunch::xp] = xp;
             particles[part][Bunch::yp] = yp;
        }
     } else {
         double reference_momentum = bunch.get_reference_particle().get_momentum();
         double m = bunch.get_mass();
         double reference_cdt = get_reference_cdt(length, k,
                                                  bunch.get_reference_particle());
         double substep_reference_cdt = reference_cdt/steps/drifts_per_step;
         double step_length = length/steps;
         double step_strength[2] = { k[0]*step_length, k[1]*step_length };

         for (int part = 0; part < local_num; ++part) {
             double x(particles[part][Bunch::x]);
             double xp(particles[part][Bunch::xp]);
             double y(particles[part][Bunch::y]);
             double yp(particles[part][Bunch::yp]);
             double cdt(particles[part][Bunch::cdt]);
             double dpop(particles[part][Bunch::dpop]);

             FF_algorithm::yoshida4<double, FF_sextupole::thin_sextupole_unit<double>, 1 >
                     ( x, xp, y, yp, cdt, dpop,
                       reference_momentum, m,
                       substep_reference_cdt,
                       step_length, step_strength, steps );

             particles[part][Bunch::x] = x;
             particles[part][Bunch::xp] = xp;
             particles[part][Bunch::y] = y;
             particles[part][Bunch::yp] = yp;
             particles[part][Bunch::cdt] = cdt;
        }
        bunch.get_reference_particle().increment_trajectory(length);
    }
}

template<class Archive>
    void
    FF_sextupole::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_sextupole::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_sextupole::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_sextupole::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_sextupole::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_sextupole::~FF_sextupole()
{

}

