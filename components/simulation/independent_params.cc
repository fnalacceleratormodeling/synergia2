#include "independent_params.h"

Independent_params::Independent_params(int map_order) :
    map_order(map_order)
{
    Operation_extractor_sptr mixed_chef_operation_extractor(
            new Chef_mixed_operation_extractor());

    extractor_map["default"] = mixed_chef_operation_extractor;
    extractor_map["mixed_chef"] = mixed_chef_operation_extractor;
    extractor_map["chef_propagate"] = Operation_extractor_sptr(
            new Chef_propagate_operation_extractor());
    extractor_map["chef_map"] = Operation_extractor_sptr(
            new Chef_map_operation_extractor());
}

int
Independent_params::get_map_order() const
{
    return map_order;
}

void
Independent_params::set_extractor(std::string const& name,
        Operation_extractor_sptr extractor)
{
    extractor_map[name] = extractor;
}

Operation_extractor_sptr
Independent_params::get_extractor(std::string const& name)
{
    return extractor_map[name];
}

Independent_params::~Independent_params()
{

}
