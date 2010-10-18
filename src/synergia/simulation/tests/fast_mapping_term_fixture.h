#ifndef FAST_MAPPING_TERM_FIXTURE_H_
#define FAST_MAPPING_TERM_FIXTURE_H_

#include "components/simulation/fast_mapping.h"

const int order = 2;
const double coeff = 3.1415;

struct Fast_mapping_term_fixture
{
    Fast_mapping_term_fixture() :
        fast_mapping_term(order)
    {
        indices = new int[order + 1];
        fast_mapping_term.coeff() = coeff;
        for (int i = 0; i < order + 1; ++i) {
            if (i <= 2) {
                indices[i] = 2 * i;
            } else {
                indices[i] = 0;
            }
            fast_mapping_term.index(i) = indices[i];
        }
    }
    ~Fast_mapping_term_fixture()
    {
        BOOST_TEST_MESSAGE("teardown Fast_mapping_term fixture");
        delete[] indices;
    }
    Fast_mapping_term fast_mapping_term;
    int *indices;
};

#endif /* FAST_MAPPING_TERM_FIXTURE_H_ */
