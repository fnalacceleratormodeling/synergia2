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
    virtual ~Distribution() = default;

    /// Get the next number in the sequence (between 0 and 1).
    virtual double get() = 0;

    /// Fill a one-dimensional array uniformly between min and max.
    virtual void fill_uniform(karray1d array, double min, double max) = 0;

    /// Fill a one-dimensional array with Gaussian distribution of
    /// zero mean and unit standard deviation.
    virtual void fill_unit_gaussian(karray1d array) = 0;

    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    virtual void fill_unit_disk(karray1d x_array, karray1d y_array) = 0;
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

    enum Generator { ranlxd2, mt19937 };

    /// Construct a Random_distribution.
    /// @param seed The random number seed. If seed == 0, the seed is
    /// obtained from Random_distribution::get_default_seed().
    /// @param comm Distribute the seed across the processors in this
    /// communicator.
    /// @param generator The underlying random number generator to be used.
    Random_distribution( unsigned long int seed, 
                         Commxx const & comm,
                         Generator generator = ranlxd2 );

    virtual ~Random_distribution();

    /// Generate a random seed. Attempt to read from device if present.
    /// Otherwise, use the system clock.
    /// @param device Read from pathname device.
    static unsigned long int get_default_seed(const char * device = "/dev/urandom");

    /// Set the random number generator seed.
    /// @param seed The seed.
    void set_seed(unsigned long int seed);

    /// Get the seed used to start the random number generator.
    unsigned long int get_original_seed() const;

    /// Get the next random number between 0 and 1.
    double get() override;

    /// Fill a one-dimensional array uniformly between min and max.
    void fill_uniform(karray1d array, double min, double max) override;
    using Distribution::fill_uniform;

    /// Fill a one-dimensional array with Gaussian distribution of
    /// zero mean and unit standard deviation.
    void fill_unit_gaussian(karray1d array) override;
    using Distribution::fill_unit_gaussian;


    /// Fill two one-dimensional arrays such that (x,y) are distributed
    /// uniformly in the unit disk.
    void fill_unit_disk(karray1d x_array, karray1d y_array) override;
    using Distribution::fill_unit_disk;
};

#endif /* DISTRIBUTION_H_ */
