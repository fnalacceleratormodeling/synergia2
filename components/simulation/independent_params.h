#ifndef INDEPENDENT_PARAMS_H_
#define INDEPENDENT_PARAMS_H_

#include "components/simulation/operation_extractor.h"
#include <string>

class Independent_params
{
private:
    int map_order;
    Operation_extractor_map extractor_map;

public:
    Independent_params(int map_order);
    int
    get_map_order() const;
    void
    set_extractor(std::string const& name, Operation_extractor_sptr extractor);
    Operation_extractor_sptr
    get_extractor(std::string const& name);
    ~Independent_params();
};

#endif /* INDEPENDENT_PARAMS_H_ */
