#include "independent_operation.h"

Independent_operation::Independent_operation(std::string const& type) :
    type(type)
{
}

Independent_operation::Independent_operation()
{
}

std::string
Independent_operation::get_type() const
{
    return type;
}

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
Fast_mapping_operation::apply(Bunch & bunch)
{
    // bunch.convert_to_state(Bunch::fixed_t);
     bunch.convert_to_state(Bunch::fixed_z_lab);
     mapping.apply(bunch);
}

Fast_mapping_operation::~Fast_mapping_operation()
{
}

Chef_propagate_operation::Chef_propagate_operation(
        Chef_lattice_section_sptr chef_lattice_section_sptr) :
    Independent_operation(chef_propagate_type_name), chef_propagator(
            chef_lattice_section_sptr)
{
}

void
Chef_propagate_operation::apply(Bunch & bunch)
{
    bunch.convert_to_state(Bunch::fixed_z_lab); 
    chef_propagator.apply(bunch);   
}

Chef_propagate_operation::~Chef_propagate_operation()
{
}
