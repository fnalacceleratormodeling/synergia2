#include "ff_sextupole.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

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
        FF_algorithm::yoshida<double, FF_algorithm::thin_sextupole_unit<double>, 4, 1 >
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
    throw std::runtime_error("Propagate JetParticle through a sextupole element is yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k2");
    k[1] = slice.get_lattice_element().get_double_attribute("k2s", 0.0);

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
        FF_algorithm::thin_sextupole_unit(x, xp, y, yp, kl);
    } else {
        FF_algorithm::yoshida<TJet<double>, FF_algorithm::thin_sextupole_unit<TJet<double> >, 4, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  step_reference_cdt,
                  step_length, step_strength, steps );
    }
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
#endif
}

void FF_sextupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    // strength
    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k2");
    k[1] = slice.get_lattice_element().get_double_attribute("k2s", 0.0);

    // tilting
    double tilt = slice.get_lattice_element().get_double_attribute("tilt", 0.0);
    if (tilt != 0.0)
    {
        std::complex<double> ck2(k[0], +k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -3.0*tilt));
        k[0] = ck2.real();
        k[1] = ck2.imag();
    }

    // scaling
    Reference_particle const & ref_l = get_ref_particle_from_slice(slice);
    Reference_particle const & ref_b = bunch.get_reference_particle();

    double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
    double brho_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

    double scale = brho_l / brho_b;

    k[0] *= scale;
    k[1] *= scale;

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();


    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

    if (length == 0.0) 
    {
        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize)
        {
            GSVector  x( &xa[part]);
            GSVector xp(&xpa[part]);
            GSVector  y( &ya[part]);
            GSVector yp(&ypa[part]);

            FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

            xp.store(&xpa[part]);
            yp.store(&ypa[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double  x( xa[part]);
            double xp(xpa[part]);
            double  y( ya[part]);
            double yp(ypa[part]);

            FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

            xpa[part] = xp;
            ypa[part] = yp;
        }
    } 
    else 
    {
        double reference_momentum = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k,
                                                 bunch.get_reference_particle());
        double step_reference_cdt = reference_cdt/steps;
        double step_length = length/steps;
        double step_strength[2] = { k[0]*step_length, k[1]*step_length };

        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize)
        {
            GSVector    x(   &xa[part]);
            GSVector   xp(  &xpa[part]);
            GSVector    y(   &ya[part]);
            GSVector   yp(  &ypa[part]);
            GSVector  cdt( &cdta[part]);
            GSVector dpop(&dpopa[part]);

            FF_algorithm::yoshida<GSVector, FF_algorithm::thin_sextupole_unit<GSVector>, 2, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );

               x.store(&xa[part]);
              xp.store(&xpa[part]);
               y.store(&ya[part]);
              yp.store(&ypa[part]);
             cdt.store(&cdta[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double    x(   xa[part]);
            double   xp(  xpa[part]);
            double    y(   ya[part]);
            double   yp(  ypa[part]);
            double  cdt( cdta[part]);
            double dpop(dpopa[part]);

            FF_algorithm::yoshida<double, FF_algorithm::thin_sextupole_unit<double>, 2, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );

               xa[part] = x;
              xpa[part] = xp;
               ya[part] = y;
              ypa[part] = yp;
             cdta[part] = cdt;
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

