#ifndef DENSE_MAPPING_H_
#define DENSE_MAPPING_H_

#include "synergia/simulation/fast_mapping.h"
#include "synergia/utils/multi_array_typedefs.h"

class Dense_mapping
{
private:
    MArray1d constant;
    MArray2d linear;
public:
    Dense_mapping(Fast_mapping const& fast_mapping);
    MArray1d_ref
    get_constant_term();
    MArray2d_ref
    get_linear_term();
    ~Dense_mapping();
};

#endif /* DENSE_MAPPING_H_ */
