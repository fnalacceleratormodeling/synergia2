#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include "utils/multi_array_typedefs.h"
#include "utils/commxx.h"
#include <gsl/gsl_rng.h>

class Distribution
{
public:
    virtual void
    fill_uniform(MArray1d_ref array) = 0;
    virtual void
    fill_unit_gaussian(MArray1d_ref array) = 0;
    virtual void
    fill_unit_disk(MArray1d_ref x_array, MArray1d_ref y_array) = 0;
    virtual
    ~Distribution()
    {
    }
    ;
};

class Random_distribution : public Distribution
{
private:
    gsl_rng *rng;
    const gsl_rng_type * rng_type;
    int rank;
    unsigned long int original_seed;
    unsigned long int
    get_default_seed();
public:
    enum Generator
    {
        ranlxd2, mt19937
    };
    Random_distribution(unsigned long int seed, Commxx const & comm,
            Generator generator = ranlxd2);
    void
    set_seed(unsigned long int seed);
    unsigned long int
    get_original_seed() const;
    virtual void
    fill_uniform(MArray1d_ref array, double min, double max);
    virtual void
    fill_unit_gaussian(MArray1d_ref array);
    virtual void
    fill_unit_disk(MArray1d_ref x_array, MArray1d_ref y_array);
    virtual
    ~Random_distribution();
};

#endif /* DISTRIBUTION_H_ */
