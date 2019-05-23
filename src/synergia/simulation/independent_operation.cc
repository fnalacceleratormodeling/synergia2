#include "independent_operation.h"
#include "synergia/libFF/ff_element_map.h"
#include "synergia/utils/simple_timer.h"


LibFF_operation::LibFF_operation(Lattice_element_slices const& slices) 
    : Independent_operation("LibFF")
    , element_slices(slices.size())
{
    populate_element_slices(slices);
}

void
LibFF_operation::apply(Bunch & bunch, Logger & logger)
{
    bunch.convert_to_state(Bunch::fixed_z_lab);

    for (auto const & es : element_slices)
        es.first.apply(es.second, bunch);
}




#if 0
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
    double t = simple_timer_current();
    bunch.convert_to_state(Bunch::fixed_z_lab);
    t = simple_timer_show(t, "fast_mapping_operation_apply-convert_to_state");
    if (verbosity > 4) {
        logger << "Fast_mapping_operation: length = " << mapping.get_length()
                << ", order = " << mapping.get_order() << std::endl;
    }
    mapping.apply(bunch);
    t = simple_timer_show(t, "fast_mapping_operation_apply-mapping_apply");
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
        Chef_lattice_section const& chef_lattice_section) :
                Independent_operation(chef_propagate_type_name),
                chef_propagator(
                        Chef_lattice_section_sptr(
                                new Chef_lattice_section(chef_lattice_section)))
{
}

Chef_propagate_operation::Chef_propagate_operation()
{
}

void
Chef_propagate_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t = simple_timer_current();
    bunch.convert_to_state(Bunch::fixed_z_lab);
    t = simple_timer_show(t, "chef_propagate_operation_apply-convert_to_state");
    chef_propagator.apply(bunch, verbosity, logger);
    t = simple_timer_show(t, "chef_propagate_operation_apply-chef_propagator_apply");
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

LibFF_operation::LibFF_operation(
        Lattice_element_slices const& lattice_element_slices) :
            Independent_operation(libFF_type_name),
            lattice_element_slices(lattice_element_slices)
{
}

LibFF_operation::LibFF_operation()
{
}

void
LibFF_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t = simple_timer_current();
    bunch.convert_to_state(Bunch::fixed_z_lab);
    t = simple_timer_show(t, "LibFF_operation_apply-convert_to_state");
    for(Lattice_element_slices::iterator it = lattice_element_slices.begin();
        it != lattice_element_slices.end(); ++it) {
        std::string const& type((*it)->get_lattice_element().get_type());
        the_big_giant_global_ff_element_map.get_element_type(type)->apply(**it, bunch);
    }
    t = simple_timer_show(t, "LibFF_operation_apply-element_apply");
}

template<class Archive>
    void
    LibFF_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Independent_operation);
        ar & BOOST_SERIALIZATION_NVP(lattice_element_slices);
    }

template
void
LibFF_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
LibFF_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
LibFF_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
LibFF_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

LibFF_operation::~LibFF_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(LibFF_operation);
#endif
