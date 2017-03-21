#include "ff_quadrupole.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

double FF_quadrupole::get_reference_cdt(double length, double * k,
                                        Reference_particle const&
                                        reference_particle) {
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
        FF_algorithm::yoshida4<double, FF_algorithm::thin_quadrupole_unit<double>, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  0.0,
                  step_length, step_strength, steps );
        reference_cdt = cdt - cdt_orig;
    }
    return reference_cdt;
}

void
FF_quadrupole::apply(Lattice_element_slice const& slice,
                     JetParticle& jet_particle)
{
    // length
    double length = slice.get_right() - slice.get_left();

    // strength
    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k[1] = slice.get_lattice_element().get_double_attribute("k1s", 0.0);

    // offsets
    const double xoff =
        slice.get_lattice_element().get_double_attribute("hoffset", 0.0);
    const double yoff =
        slice.get_lattice_element().get_double_attribute("voffset", 0.0);

    const double vxoff(xoff);
    const double vyoff(yoff);

    // tilting
    double tilt = slice.get_lattice_element().get_double_attribute("tilt", 0.0);
    if (tilt != 0.0) {
        std::complex<double> ck2(k[0], -k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -2.0 * tilt));
        k[0] = ck2.real();
        k[1] = ck2.imag();
    }

    // scaling
    Reference_particle const& ref_l = get_ref_particle_from_slice(slice);

    double brho_l = ref_l.get_momentum() / ref_l.get_charge(); // GV/c
    double brho_b = jet_particle.Momentum().standardPart() *
                    (1.0 + jet_particle.get_npz().standardPart()) /
                    ref_l.get_charge(); // GV/c

    double scale = brho_l / brho_b;

    k[0] *= scale;
    k[1] *= scale;

    if (length == 0.0) {
        Jet& x(jet_particle.State()[Particle::xIndex]);
        Jet& xp(jet_particle.State()[Particle::npxIndex]);
        Jet& y(jet_particle.State()[Particle::yIndex]);
        Jet& yp(jet_particle.State()[Particle::npyIndex]);

        x -= vxoff;
        y -= vyoff;

        FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k);

        x += vxoff;
        y += vyoff;
    } else {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute(
            "yoshida_steps", 4.0);

        // params
        double reference_momentum = jet_particle.Momentum().standardPart();
        double m = jet_particle.Mass();
        double reference_cdt = get_reference_cdt(length, k, ref_l);
        double step_reference_cdt = reference_cdt / steps;
        double step_length = length / steps;
        double step_strength[2] = { k[0] * step_length, k[1] * step_length };

        Jet& x(jet_particle.State()[Particle::xIndex]);
        Jet& xp(jet_particle.State()[Particle::npxIndex]);
        Jet& y(jet_particle.State()[Particle::yIndex]);
        Jet& yp(jet_particle.State()[Particle::npyIndex]);
        Jet& cdt(jet_particle.State()[Particle::cdtIndex]);
        Jet& dpop(jet_particle.State()[Particle::ndpIndex]);

        x -= vxoff;
        y -= vyoff;

#if 1
        FF_algorithm::yoshida6<Jet, FF_algorithm::thin_quadrupole_unit<Jet>, 1>(
            x, xp, y, yp, cdt, dpop, reference_momentum, m, step_reference_cdt,
            step_length, step_strength, steps);
#endif

        x += vxoff;
        y += vyoff;
    }
}

void FF_quadrupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    // length
    double length = slice.get_right() - slice.get_left();

    // strength
    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k1", 0.0);
    k[1] = slice.get_lattice_element().get_double_attribute("k1s", 0.0);

    // offsets
    const double xoff = slice.get_lattice_element().get_double_attribute("hoffset", 0.0);
    const double yoff = slice.get_lattice_element().get_double_attribute("voffset", 0.0);

    const GSVector vxoff(xoff);
    const GSVector vyoff(yoff);

    // tilting
    double tilt = slice.get_lattice_element().get_double_attribute("tilt", 0.0);
    if (tilt != 0.0)
    {
        std::complex<double> ck2(k[0], -k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -2.0*tilt));
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

            x -= vxoff;
            y -= vyoff;

            FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k);

            xp.store(&xpa[part]);
            yp.store(&ypa[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double  x( xa[part]);
            double xp(xpa[part]);
            double  y( ya[part]);
            double yp(ypa[part]);

            x -= xoff;
            y -= yoff;

            FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k);

            xpa[part] = xp;
            ypa[part] = yp;
        }
    } 
    else 
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // params
        double reference_momentum = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k, bunch.get_reference_particle());
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

            x -= vxoff;
            y -= vyoff;

#if 0
            FF_algorithm::yoshida<GSVector, FF_algorithm::thin_quadrupole_unit<GSVector>, 6/*order*/, 1/*components*/ >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );
#endif

#if 1
            FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_quadrupole_unit<GSVector>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );
#endif

#if 0
            FF_algorithm::yoshida4<GSVector, FF_algorithm::thin_quadrupole_unit<GSVector>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );
#endif

#if 0
            FF_algorithm::quadrupole_chef<GSVector>
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      reference_cdt, length, k, 40 /* kicks */ );
#endif

            x += vxoff;
            y += vyoff;

               x.store(&xa[part]);
              xp.store(&xpa[part]);
               y.store(&ya[part]);
              yp.store(&ypa[part]);
             cdt.store(&cdta[part]);
            dpop.store(&dpopa[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double    x(   xa[part]);
            double   xp(  xpa[part]);
            double    y(   ya[part]);
            double   yp(  ypa[part]);
            double  cdt( cdta[part]);
            double dpop(dpopa[part]);

            x -= xoff;
            y -= yoff;

            FF_algorithm::yoshida6<double, FF_algorithm::thin_quadrupole_unit<double>, 1 >
                    ( x, xp, y, yp, cdt, dpop,
                      reference_momentum, m,
                      step_reference_cdt,
                      step_length, step_strength, steps );

            x += xoff;
            y += yoff;

               xa[part] = x;
              xpa[part] = xp;
               ya[part] = y;
              ypa[part] = yp;
             cdta[part] = cdt;
            dpopa[part] = dpop;
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
}

template<class Archive>
    void
    FF_quadrupole::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_quadrupole::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_quadrupole::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_quadrupole::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_quadrupole::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_quadrupole::~FF_quadrupole()
{

}

