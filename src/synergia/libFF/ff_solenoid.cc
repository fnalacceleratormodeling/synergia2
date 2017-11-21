#include "ff_solenoid.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

FF_solenoid::FF_solenoid()
{
}

double FF_solenoid::get_reference_cdt_drift(double length, Reference_particle & reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    double reference_momentum = reference_particle.get_momentum()*(1+dpop);
    double m = reference_particle.get_mass();

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m, 0.0);

    // propagate and update the bunch design reference particle state
    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}

double FF_solenoid::get_reference_cdt_solenoid(double length, Reference_particle & reference_particle,
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


void FF_solenoid::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a solenoid element is yet to be implemented");
}

void FF_solenoid::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const int  local_num = bunch.get_local_num();
    const double    mass = bunch.get_mass();

    // in edge, out edge
    const bool has_in_edge  = slice.has_left_edge();
    const bool has_out_edge = slice.has_right_edge();

    // ks here is the field
    double ks = slice.get_lattice_element().get_double_attribute("ks");

    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();
    const double   ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

    if (fabs(ks) < 1e-12)
    {
        // reference cdt
        const double ref_cdt = get_reference_cdt_drift(length, ref_l);

        // this is a drift
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
    }
    else
    {
        // scaling
        double brho_l = ref_l.get_momentum() * (1.0 + ref_l.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c
        double brho_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

        double scale = brho_l / brho_b;

                ks = ks * scale;
        double ksl = ks * length;
        double kse = ks * 0.5;

        // reference cdt
        const double ref_cdt = get_reference_cdt_solenoid(length, ref_l, 
                has_in_edge, has_out_edge, ks, kse, ksl);

        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part)
        {
            // in edge
            if (has_in_edge)
            {
                FF_algorithm::solenoid_in_edge_kick(xa[part], xpa[part], ya[part], ypa[part], kse);
            }

            // body
            FF_algorithm::solenoid_unit(
                    xa[part], xpa[part], ya[part], ypa[part], cdta[part], dpopa[part],
                    ksl, ks, length, ref_p, mass, ref_cdt );

            // out edge
            if (has_out_edge)
            {
                FF_algorithm::solenoid_out_edge_kick(xa[part], xpa[part], ya[part], ypa[part], kse);
            }
        }

    }

    bunch.get_reference_particle().increment_trajectory(length);
}

template<class Archive>
    void
    FF_solenoid::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_solenoid::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_solenoid::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_solenoid::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_solenoid::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_solenoid::~FF_solenoid()
{

}
