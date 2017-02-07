#include "ff_constfoc.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"
#include "synergia/foundation/math_constants.h"

FF_constfoc::FF_constfoc()
{
}

double FF_constfoc::get_reference_cdt(
        double length, double csl, double snl, double BL, double iBL,
        Reference_particle & reference_particle)
{
    double cdt(reference_particle.get_state()[Bunch::cdt]);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    double cdt_orig = cdt;
    FF_algorithm::constfoc_unit(cdt, dpop, csl, snl, BL, iBL);

    return cdt - cdt_orig;
}

void FF_constfoc::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("FF_constfoc for Jet Particles are not implemented yet");
}

void FF_constfoc::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const int  local_num = bunch.get_local_num();

    const double BH = slice.get_lattice_element().get_double_attribute("betah");
    const double BV = slice.get_lattice_element().get_double_attribute("betav");
    const double BL = slice.get_lattice_element().get_double_attribute("betal");
    const double PL = slice.get_lattice_element().get_double_attribute("nul") * 2.0 * mconstants::pi;

    const double iBH = 1.0 / BH;
    const double iBV = 1.0 / BV;
    const double iBL = 1.0 / BL;

    const double csh = cos(length * iBH);
    const double snh = sin(length * iBH);

    const double csv = cos(length * iBV);
    const double snv = sin(length * iBV);

    const double csl = cos(PL);
    const double snl = sin(PL);

    const double ref_cdt = get_reference_cdt(length, csl, snl, BL, iBL, bunch.get_reference_particle());

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

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

        FF_algorithm::constfoc_unit(x, xp, csh, snh, BH, iBH);
        FF_algorithm::constfoc_unit(y, yp, csv, snv, BV, iBV);
        FF_algorithm::constfoc_unit(cdt, dpop, csl, snl, BL, iBL);

        x.store(&xa[part]);
        xp.store(&xpa[part]);
        y.store(&ya[part]);
        yp.store(&ypa[part]);
        cdt.store(&cdta[part]);
        dpop.store(&dpopa[part]);

        for(int i=0; i<gsvsize; ++i)
        {
            cdta[part+i] -= ref_cdt;
        }
    }

    for (int part = block_last; part < local_num; ++part) 
    {
        double x(xa[part]);
        double xp(xpa[part]);
        double y(ya[part]);
        double yp(ypa[part]);
        double cdt(cdta[part]);
        double dpop(dpopa[part]);

        FF_algorithm::constfoc_unit(x, xp, csh, snh, BH, iBH);
        FF_algorithm::constfoc_unit(y, yp, csv, snv, BV, iBV);
        FF_algorithm::constfoc_unit(cdt, dpop, csl, snl, BL, iBL);

        xa[part] = x;
        xpa[part] = xp;
        ya[part] = y;
        ypa[part] = yp;
        cdta[part] = cdt - ref_cdt;
        dpopa[part] = dpop;
    }

    bunch.get_reference_particle().increment_trajectory(length);
}

template<class Archive>
    void
    FF_constfoc::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_constfoc::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_constfoc::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_constfoc::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_constfoc::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_constfoc::~FF_constfoc()
{

}
