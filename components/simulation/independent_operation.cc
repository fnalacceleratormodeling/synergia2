#include "independent_operation.h"

Fast_mapping_operation::Fast_mapping_operation(Fast_mapping const& mapping) :
    mapping(mapping)
{
}

void
Fast_mapping_operation::apply(Bunch & bunch)
{
    mapping.apply(bunch);
}

Fast_mapping_operation::~Fast_mapping_operation()
{
}

Chef_propagate_operation::Chef_propagate_operation(
        Chef_elements const& chef_elements) :
    chef_propagator(chef_elements)
{
}

void
Chef_propagate_operation::apply(Bunch & bunch)
{
    chef_propagator.apply(bunch);
}

Chef_propagate_operation::~Chef_propagate_operation()
{
}
