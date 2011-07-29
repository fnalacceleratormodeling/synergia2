#include "independent_operation.h"

Independent_operation::Independent_operation(std::string const& type) :
    type(type)
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

void
Fast_mapping_operation::apply(Bunch & bunch)
{ 
     bunch.convert_to_state(Bunch::fixed_t);
     bunch.convert_to_state(Bunch::fixed_z);     
     mapping.apply(bunch);
}

Fast_mapping_operation::~Fast_mapping_operation()
{
}

Chef_propagate_operation::Chef_propagate_operation(
        Chef_elements const& chef_elements) :
    Independent_operation(chef_propagate_type_name), chef_propagator(
            chef_elements)
{
}

void
Chef_propagate_operation::apply(Bunch & bunch)
{
    bunch.convert_to_state(Bunch::fixed_z); 
    chef_propagator.apply(bunch);   
}

Chef_propagate_operation::~Chef_propagate_operation()
{
}
