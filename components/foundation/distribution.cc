#include "distribution.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <fstream>
#include <stdexcept>

void
Distribution::fill_uniform(MArray1d_ref array, double min, double max)
{
    fill_uniform(array[boost::indices[range()]], min, max);
}

void
Distribution::fill_unit_gaussian(MArray1d_ref array)
{
    fill_unit_gaussian(array[boost::indices[range()]]);
}

void
Distribution::fill_unit_disk(MArray1d_ref x_array, MArray1d_ref y_array)
{
    fill_unit_disk(x_array[boost::indices[range()]],
            y_array[boost::indices[range()]]);
}

unsigned long int
Random_distribution::get_default_seed(const char * device)
{
    unsigned long int seed;
    std::ifstream devrandom(device, std::ios::binary);
    if (devrandom) {
        devrandom.read((char *) &seed, sizeof(unsigned long int));
    } else {
        seed = std::time(0);
    }
    return seed;
}

Random_distribution::Random_distribution(unsigned long int seed,
        Commxx const & comm, Generator generator)
{
    gsl_rng_env_setup();
    if (generator == ranlxd2) {
        rng_type = gsl_rng_ranlxd2;
    } else if (generator == mt19937) {
        rng_type = gsl_rng_mt19937;
    }
    rng = gsl_rng_alloc(rng_type);
    rank = comm.get_rank();
    set_seed(seed);
}

void
Random_distribution::set_seed(unsigned long int seed)
{
    if (seed == 0) {
        original_seed = get_default_seed();
    } else {
        original_seed = seed;
    }
    unsigned long int distributed_seed;
    distributed_seed = (1000 + 5 * (rank + original_seed)) * ((rank
            + original_seed) + 7) - 1;

    gsl_rng_set(rng, distributed_seed);
}

unsigned long int
Random_distribution::get_original_seed() const
{
    return original_seed;
}

double
Random_distribution::get()
{
    return gsl_rng_uniform(rng);
}

void
Random_distribution::fill_uniform(MArray1d_view array, double min, double max)
{
    for (MArray1d::iterator it = array.begin(); it != array.end(); ++it) {
        *it = gsl_ran_flat(rng, min, max);
    }
}

void
Random_distribution::fill_unit_gaussian(MArray1d_view array)
{
    for (MArray1d::iterator it = array.begin(); it != array.end(); ++it) {
        *it = gsl_ran_ugaussian_ratio_method(rng);
    }

}

void
Random_distribution::fill_unit_disk(MArray1d_view x_array,
        MArray1d_view y_array)
{
    if (x_array.shape()[0] != y_array.shape()[0]) {
        throw std::runtime_error(
                "Random_distribution::fill_unit_disk: x_array and y_array must have the same length\n");
    }
    for (unsigned int n = 0; n < x_array.shape()[0]; ++n) {
        double x = 1.0, y = 1.0;
        while (x * x + y * y > 1.0) {
            x = 2.0 * gsl_rng_uniform(rng) - 1.0;
            y = 2.0 * gsl_rng_uniform(rng) - 1.0;
        }
        x_array[n] = x;
        y_array[n] = y;
    }
}

Random_distribution::~Random_distribution()
{
    gsl_rng_free(rng);
}
