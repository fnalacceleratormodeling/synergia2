#include "independent_operation.h"
#include "synergia/utils/simple_timer.h"

Independent_operation::Independent_operation(std::string const& type) :
    type(type)
{
}

Independent_operation::Independent_operation()
{
}

std::string const&
Independent_operation::get_type() const
{
    return type;
}

template<class Archive>
    void
    Independent_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(type);
    }

template
void
Independent_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Independent_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Independent_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Independent_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Independent_operation::~Independent_operation()
{
}

Fast_mapping_operation::Fast_mapping_operation(Fast_mapping const& mapping) :
    Independent_operation(fast_mapping_type_name), mapping(mapping)
{
}

Fast_mapping_operation::Fast_mapping_operation()
{
}

void
Fast_mapping_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    bunch.convert_to_state(Bunch::fixed_z_lab);
    if (verbosity > 4) {
        logger << "Fast_mapping_operation: length = " << mapping.get_length()
                << ", order = " << mapping.get_order() << std::endl;
    }
    mapping.apply(bunch);
}

Fast_mapping const&
Fast_mapping_operation::get_fast_mapping() const
{
    return mapping;
}

template<class Archive>
    void
    Fast_mapping_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Independent_operation);
        ar & BOOST_SERIALIZATION_NVP(mapping);
    }

template
void
Fast_mapping_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Fast_mapping_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Fast_mapping_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Fast_mapping_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Fast_mapping_operation::~Fast_mapping_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Fast_mapping_operation);

Chef_propagate_operation::Chef_propagate_operation(
        Chef_lattice_section_sptr chef_lattice_section_sptr) :
    Independent_operation(chef_propagate_type_name), chef_propagator(
            chef_lattice_section_sptr)
{
}

Chef_propagate_operation::Chef_propagate_operation()
{
}

void
Chef_propagate_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    bunch.convert_to_state(Bunch::fixed_z_lab);
    chef_propagator.apply(bunch, verbosity, logger);
}

template<class Archive>
    void
    Chef_propagate_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Independent_operation);
        ar & BOOST_SERIALIZATION_NVP(chef_propagator);
    }

template
void
Chef_propagate_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Chef_propagate_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Chef_propagate_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Chef_propagate_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Chef_propagate_operation::~Chef_propagate_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Chef_propagate_operation);
