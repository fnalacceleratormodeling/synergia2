#include "ff_mcmillan.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

double FF_mcmillan::get_reference_cdt(double length, double j0, double beta_e, double radius, Reference_particle &reference_particle, bool simple)
{
    double pref = reference_particle.get_momentum();
    double m = reference_particle.get_mass();
    double beta_p = reference_particle.get_beta();

    double x   (reference_particle.get_state()[Bunch::x]);
    double xp  (reference_particle.get_state()[Bunch::xp]);
    double y   (reference_particle.get_state()[Bunch::y]);
    double yp  (reference_particle.get_state()[Bunch::yp]);
    double cdt (0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    if ( close_to_zero(length) )
    {
        // for 0 length, j0 is the integrated current density*length
        FF_algorithm::thin_mcmillan_unit(x, xp, y, yp, dpop, j0, beta_p, beta_e, radius);
    }
    else
    {
        if (simple)
        {
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, 0.0);

            // for >0 length, j0 is the current density, so must be multiplied by length for integrated strength
            FF_algorithm::thin_mcmillan_unit(x, xp, y, yp, dpop, j0*length, beta_p, beta_e, radius);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, 0.0);
        }
        else
        {
            // steps comes from base class, set in apply method
            double step_length = length / steps;

            // for >0 length, j0 is the current density
            double step_strength = j0*step_length;

            // propagate
            FF_algorithm::yoshida6<double, FF_algorithm::thin_mcmillan_unit<double>, 1>
                ( x, xp, y, yp, cdt, dpop,
                  pref, m, 0.0,
                  step_length, step_strength, steps );
        }
    }

    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}

void FF_mcmillan::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF::mcmillan for JetParticles has yet to be implemented");

}


void FF_mcmillan::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    double   l = slice.get_lattice_element().get_double_attribute("l", 0.0);
    double j0 = slice.get_lattice_element().get_double_attribute("j0");
    double beta_e = slice.get_lattice_element().get_double_attribute("ebeta");
    double radius = slice.get_lattice_element().get_double_attribute("radius");

    Reference_particle       & ref_lattice = bunch.get_design_reference_particle();
    Reference_particle const & ref_bunch   = bunch.get_reference_particle();

    double plattice = ref_lattice.get_momentum();
    double pbunch = ref_bunch.get_momentum();

    double beta_p = ref_bunch.get_beta();
    double scale = ( ref_bunch.get_charge() / ref_lattice.get_charge() );
    j0 *= scale; // account for possible multiply charged beam particles

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    const int gsvsize = GSVector::size();

    double * RESTRICT xa;
    double * RESTRICT xpa;
    double * RESTRICT ya;
    double * RESTRICT ypa;
    double * RESTRICT cdta;
    double * RESTRICT dpopa;

    if ( close_to_zero(length) )
    {
        // only to update the reference particle
        // the reference time is calculated with the design reference particle which
        // is at plattice
        double reference_cdt = get_reference_cdt(l, j0, beta_e, radius. ref_bunch, false);

        // real particles
        {
            bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

            const int num_blocks = local_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector x(&xa[part]);
                GSVector y(&ya[part]);
                GSVector xp(&xpa[part]);
                GSVector yp(&ypa[part]);
                GSVector dpop(&dpopa[part]);

                FF_algorithm::thin_mcmillan_unit(x, xp, y, yp, dpop, j0, beta_p, beta_e, radius);

                x.store(&xa[part]);
                y.store(&ya[part]);
                xp.store(&xpa[part]);
                yp.store(&ypa[part]);
                dpopa.store(&dpop[part]);
            }

            for (int part = block_last; part < local_num; ++part)
            {
                double x(xa[part]);
                double y(ya[part]);
                double xp(xpa[part]);
                double yp(ypa[part]);
                double dpop(dpopa[part]);

                FF_algorithm::thin_mcmillan_unit(x. xp, y, yp, dpop, j0, beta_p, beta_e, radius);

                xa[part] = x;
                ya[part] = y;
                xpa[part] = xp;
                ypa[part] = yp;
                dpopa[part] = dpop;
            }
        }

        // spectator particles
        {
            bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

            const int num_blocks = local_s_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector x(&xa[part]);
                GSVector y(&ya[part]);
                GSVector xp(&xpa[part]);
                GSVector yp(&ypa[part]);
                GSVector dpop(&dpopa[part]);

                FF_algorithm::thin_mcmillan_unit(x, xp, y, yp, dpop, j0, beta_p, beta_e, radius);

                x.store(&xa[part]);
                y.store(&ya[part]);
                xp.store(&xpa[part]);
                yp.store(&ypa[part]);
                dpopa.store(&dpop[part]);
            }

            for (int part = block_last; part < local_s_num; ++part)
            {
                double x(xa[part]);
                double y(ya[part]);
                double xp(xpa[part]);
                double yp(ypa[part]);
                double dpop(dpopa[part]);

                FF_algorithm::thin_mcmillan_unit(x, xp, y, yp, dpop, j0, beta_p, beta_e, radius);

                xa[part] = x;
                ya[part] = y;
                xpa[part] = xp;
                ypa[part] = yp;
                dpopa[part] = dpop;
            }
        }
    }
    else
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // simple drift-kick-drift scheme
        double simple_d = slice.get_lattice_element().get_double_attribute("simple", 0.0);
        bool simple = fabs(simple_d) > 1e-16;

        // strength per unit length, pul = per unit length
        double j0_pul = j0 / l;

        double pref = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(l, j0, beta_e, ref_bunch, simple);

        if (simple)
        {
            // real particles
            {
                bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

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

                    // simple drift - kick - drift scheme, this is how CHEF does it
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_mcmillan_unit(xp, b_hk);
                    FF_algorithm::thin_mcmillan_unit(yp, b_vk);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                       x.store(&xa[part]);
                      xp.store(&xpa[part]);
                       y.store(&ya[part]);
                      yp.store(&ypa[part]);
                     cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_mcmillan_unit(xp, b_hk);
                    FF_algorithm::thin_mcmillan_unit(yp, b_vk);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                       xa[part] = x;
                      xpa[part] = xp;
                       ya[part] = y;
                      ypa[part] = yp;
                     cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }

            // spectator particles
            {
                bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

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

                    // simple drift - kick - drift scheme, this is how CHEF does it
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_mcmillan_unit(xp, b_hk);
                    FF_algorithm::thin_mcmillan_unit(yp, b_vk);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                       x.store(&xa[part]);
                      xp.store(&xpa[part]);
                       y.store(&ya[part]);
                      yp.store(&ypa[part]);
                     cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_s_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_mcmillan_unit(xp, b_hk);
                    FF_algorithm::thin_mcmillan_unit(yp, b_vk);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                       xa[part] = x;
                      xpa[part] = xp;
                       ya[part] = y;
                      ypa[part] = yp;
                     cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }
        }
        else
        {
            double b_hk_pul = b_hk / l;
            double b_vk_pul = b_vk / l;

            double step_reference_cdt = reference_cdt / steps;
            double step_length = length / steps;
            double step_strength[2] = { b_hk_pul * step_length, b_vk_pul * step_length };

            // real particles
            {
                bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

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

                    // yoshida kick
                    FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_mcmillan_unit<GSVector>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                       x.store(&xa[part]);
                      xp.store(&xpa[part]);
                       y.store(&ya[part]);
                      yp.store(&ypa[part]);
                     cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::yoshida6<double, FF_algorithm::thin_mcmillan_unit<double>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                       xa[part] = x;
                      xpa[part] = xp;
                       ya[part] = y;
                      ypa[part] = yp;
                     cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }

            // spectator particles
            {
                bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

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

                    // yoshida kick
                    FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_mcmillan_unit<GSVector>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                       x.store(&xa[part]);
                      xp.store(&xpa[part]);
                       y.store(&ya[part]);
                      yp.store(&ypa[part]);
                     cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_s_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::yoshida6<double, FF_algorithm::thin_mcmillan_unit>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                       xa[part] = x;
                      xpa[part] = xp;
                       ya[part] = y;
                      ypa[part] = yp;
                     cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
}

template<class Archive>
    void
    FF_mcmillan::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_mcmillan::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_mcmillan::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_mcmillan::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_mcmillan::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_mcmillan::~FF_mcmillan()
{

}

