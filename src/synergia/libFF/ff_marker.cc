#include "ff_marker.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"
#include <iomanip>

FF_marker::FF_marker()
{
}

double FF_marker::get_reference_cdt(double length,
                                   Reference_particle & reference_particle)
{
    // marker does nothing
    return 0.0;
}

void FF_marker::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    return;
}

void FF_marker::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    // marker does nothing
    return;
}

template<class Archive>
    void
    FF_marker::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_marker::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_marker::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_marker::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_marker::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_marker::~FF_marker()
{

}
