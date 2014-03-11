#ifndef DENSE_MAPPING_H_
#define DENSE_MAPPING_H_

#include "synergia/simulation/fast_mapping.h"
#include "synergia/utils/multi_array_typedefs.h"

class Dense_mapping
{
private:
    double length;
    MArray1d constant;
    MArray2d linear;
public:
    Dense_mapping();
    Dense_mapping(Fast_mapping const& fast_mapping);
    double
    get_length() const;
    MArray1d_ref
    get_constant_term() const;
    MArray2d_ref
    get_linear_term() const;
    ~Dense_mapping();
};

#endif /* DENSE_MAPPING_H_ */
