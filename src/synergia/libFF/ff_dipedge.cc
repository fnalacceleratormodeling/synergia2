#include "ff_dipedge.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

FF_dipedge::FF_dipedge()
{
}

double FF_dipedge::get_reference_cdt(double re_2_1, double re_4_3, double * te, Reference_particle & reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    // propagate
    FF_algorithm::dipedge_unit(x, xp, y, yp, re_2_1, re_4_3, te);

    // update the bunch design reference particle state
    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}

void FF_dipedge::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a dipedge element is yet to be implemented");
}

void FF_dipedge::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double    h = slice.get_lattice_element().get_double_attribute("h");
    const double edge = slice.get_lattice_element().get_double_attribute("e1");
    const double fint = slice.get_lattice_element().get_double_attribute("fint");
    const double hgap = slice.get_lattice_element().get_double_attribute("hgap");
    const double tilt = slice.get_lattice_element().get_double_attribute("tilt");
    const double corr = (h + h) * hgap * fint;

    const bool   fsec = false;
    const double  sig = 1.0;
    const double   he = 0.0;
    const double  sk1 = 0.0;

    // -----------
    // Linear part

    double tanedg =   tan(edge);
    double secedg =   1.0 / cos(edge);
    double sn     =   sin(edge);
    double psip   =   edge - corr * secedg * (1.0 + sn * sn);
    double re_2_1 =   h * tanedg;
    double re_4_3 = - h * tan(psip);

    // ----------------------
    // Second order terms

    // te[0] = te_1_1_1, te[1] = te_1_3_3
    // te[2] = te_2_1_1, te[3] = te_2_1_2, te[4] = te_2_3_3, te[5] = te_2_3_4
    // te[6] = te_3_1_3
    // te[7] = te_4_1_3, te[8] = te_4_1_4, te[9] = te_4_2_3
    double te[10] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    if (fsec) 
    {
        double hh = sig * ( h / 2.0 );

        te[0] = -hh*tanedg*tanedg;
        te[1] =  hh*secedg*secedg;

        te[2] =   (h/2.0) * he * (secedg*secedg*secedg) + sk1 * tanedg;
        te[3] = - te[0];
        te[4] =   hh * h * (tanedg*tanedg*tanedg)  -  te[2];
        te[5] =   te[0];

        te[6] = - te[0];

        te[7] = - te[2];
        te[8] =   te[0];
        te[9] = - te[1];

        if( sig > 0 ) 
        {
            te[4] += (h*secedg)*(h*secedg)*(tanedg/2.0);
        }
        else 
        {
            te[2] -= (h*secedg)*(h*secedg)*(tanedg/2.0);
            te[7] += (h*secedg)*(h*secedg)*(tanedg/2.0);
        }
    }

    // propagate the reference particle
    Reference_particle & ref = bunch.get_design_reference_particle();
    get_reference_cdt(re_2_1, re_4_3, te, ref);

    // propagate bunch particles
    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    const int local_num = bunch.get_local_num();
    const int local_s_num = bunch.get_local_spectator_num();

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

            FF_algorithm::dipedge_unit(x, xp, y, yp, re_2_1, re_4_3, te);

            x.store(&xa[part]);
            xp.store(&xpa[part]);
            y.store(&ya[part]);
            yp.store(&ypa[part]);
        }
     
        for (int part = block_last; part < local_num; ++part)
        {
            FF_algorithm::dipedge_unit(xa[part], xpa[part], ya[part], ypa[part], re_2_1, re_4_3, te);
        }
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

            FF_algorithm::dipedge_unit(x, xp, y, yp, re_2_1, re_4_3, te);

            x.store(&xa[part]);
            xp.store(&xpa[part]);
            y.store(&ya[part]);
            yp.store(&ypa[part]);
        }
     
        for (int part = block_last; part < local_s_num; ++part)
        {
            FF_algorithm::dipedge_unit(xa[part], xpa[part], ya[part], ypa[part], re_2_1, re_4_3, te);
        }
    }
}

template<class Archive>
    void
    FF_dipedge::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_dipedge::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_dipedge::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_dipedge::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_dipedge::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_dipedge::~FF_dipedge()
{

}
