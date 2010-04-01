#include "independent_operation.h"
#include "chef_propagate.h"

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
    chef_elements(chef_elements)
{
}

void
Chef_propagate_operation::apply(Bunch & bunch)
{
    chef_propagate(bunch, chef_elements);
}

Chef_propagate_operation::~Chef_propagate_operation()
{
}
