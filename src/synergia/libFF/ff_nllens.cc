#include "ff_nllens.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

FF_nllens::FF_nllens()
{
}

double FF_nllens::get_reference_cdt(double icnll, double kick, Reference_particle & reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    // propagate
    FF_algorithm::nllens_unit(x, y, xp, yp, icnll, kick);

    if (std::isnan(xp) || std::isnan(yp))
        throw std::runtime_error("the bunch reference particle hits the sigularity region when propagating through a nonlinearlens");

    // update the bunch design reference particle state
    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}

void FF_nllens::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("Propagate JetParticle through a nllens element is yet to be implemented");
}

void FF_nllens::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double knll = slice.get_lattice_element().get_double_attribute("knll");
    const double cnll = slice.get_lattice_element().get_double_attribute("cnll");

    const double icnll = 1.0 / cnll;
    const double kick = -knll * icnll;

    Reference_particle & ref = bunch.get_design_reference_particle();
    get_reference_cdt(icnll, kick, ref);

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    // bunch particles
    {
        bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);
        const int  local_num = bunch.get_local_num();

        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part)
        {
            FF_algorithm::nllens_unit(xa[part], ya[part], xpa[part], ypa[part], icnll, kick);
        }
    }

    // bunch spectator particles
    {
        bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);
        const int  local_s_num = bunch.get_local_spectator_num();

        #pragma omp parallel for
        for (int part = 0; part < local_s_num; ++part)
        {
            FF_algorithm::nllens_unit(xa[part], ya[part], xpa[part], ypa[part], icnll, kick);
        }
    }
}

template<class Archive>
    void
    FF_nllens::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_nllens::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_nllens::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_nllens::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_nllens::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_nllens::~FF_nllens()
{

}
