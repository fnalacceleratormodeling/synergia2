#include "dense_mapping.h"

Dense_mapping::Dense_mapping(Fast_mapping const& fast_mapping) :
    constant(boost::extents[6]), linear(boost::extents[6][6])
{
    for (int i = 0; i < 6; ++i) {
        constant[i] = 0.0;
        for (int j = 0; j < 6; ++j) {
            linear[i][j] = 0.0;
        }
    }
    std::vector<std::vector<Fast_mapping_terms > > const& terms =
            fast_mapping.get_terms();
    for (int i = 0; i < 6; ++i) {
        Fast_mapping_terms::const_iterator telem;
        for (telem = terms[i][0].begin(); telem != terms[i][0].end(); ++telem) {
            constant[i] += telem->coeff();
        }
        for (telem = terms[i][1].begin(); telem != terms[i][1].end(); ++telem) {
            linear[i][telem->index(0)] = telem->coeff();
        }
    }
}

MArray1d_ref
Dense_mapping::get_constant_term()
{
    return constant;
}

MArray2d_ref
Dense_mapping::get_linear_term()
{
    return linear;
}

Dense_mapping::~Dense_mapping()
{

}

