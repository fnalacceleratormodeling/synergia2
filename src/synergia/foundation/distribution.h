#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"
#include <gsl/gsl_rng.h>

/// Distribution is a virtual base class for obtaining the next number or set
/// of numbers from a sequence according to a limited set of shapes.
class Distribution
{
public:
    /// Get the next number in the sequence (between 0 and 1).
    virtual double
    get() = 0;

    /// Fill a one-dimensional array uniformly between min and max.
    virtual void
    fill_uniform(MArray1d_view array, double min, double max) = 0;

    /// Alternate form for type compatibility.
    void
    fill_uniform(MArray1d_ref array, double min, double max);

    /// Fill a one-dimensional array with Gaussian distribution of
    /// zero mean and unit standard deviation.
    virtual void
    fill_unit_gaussian(MArray1d_view array) = 0;

    /// Alternate form for type compatibility.
    void
    fill_unit_gaussian(MArray1d_ref array);

    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    virtual void
    fill_unit_disk(MArray1d_view x_array, MArray1d_view y_array) = 0;

    /// Alternate form for type compatibility.
    void
    fill_unit_disk(MArray1d_ref x_array, MArray1d_ref y_array);

    virtual
    ~Distribution()
    {
    }
    ;
};

/// Random_distribution provides a Distribution of random numbers. The random seed
/// is maintained across multiple processors. The implementation uses random numbers
/// from the GNU Scientific Library.
class Random_distribution : public Distribution
{
private:
    gsl_rng *rng;
    const gsl_rng_type * rng_type;
    int rank;
    unsigned long int original_seed;
public:
    enum Generator
    {
        ranlxd2, mt19937
    };

    /// Construct a Random_distribution.
    /// @param seed The random number seed. If seed == 0, the seed is
    /// obtained from Random_distribution::get_default_seed().
    /// @param comm Distribute the seed across the processors in this
    /// communicator.
    /// @param generator The underlying random number generator to be used.
    Random_distribution(unsigned long int seed, Commxx const & comm,
            Generator generator = ranlxd2);

    /// Generate a random seed. Attempt to read from device if present.
    /// Otherwise, use the system clock.
    /// @param device Read from pathname device.
    static unsigned long int
    get_default_seed(const char * device = "/dev/urandom");

    /// Set the random number generator seed.
    /// @param seed The seed.
    void
    set_seed(unsigned long int seed);

    /// Get the seed used to start the random number generator.
    unsigned long int
    get_original_seed() const;

    /// Get the next random number between 0 and 1.
    virtual double
    get();

    /// Fill a one-dimensional array uniformly between min and max.
    virtual void
    fill_uniform(MArray1d_view array, double min, double max);

    using Distribution::fill_uniform;

    /// Fill a one-dimensional array with Gaussian distribution of
    /// zero mean and unit standard deviation.
    virtual void
    fill_unit_gaussian(MArray1d_view array);

    using Distribution::fill_unit_gaussian;

    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    virtual void
    fill_unit_disk(MArray1d_view x_array, MArray1d_view y_array);

    using Distribution::fill_unit_disk;

    virtual
    ~Random_distribution();
};

#endif /* DISTRIBUTION_H_ */
