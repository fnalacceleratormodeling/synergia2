#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

#include "boost/multi_array.hpp"

class Reference_particle
{
private:
    double total_energy;
    boost::multi_array<double, 1 > state;
public:
    Reference_particle(double total_energy);
    Reference_particle(double total_energy, boost::const_multi_array_ref<
            double, 1 > state);

    void
    set_total_energy(double total_energy);
    void
    set_state(boost::const_multi_array_ref<double, 1 > state);

    double
    get_total_energy();
    boost::multi_array_ref<double, 1 >
    get_state();
};

#endif /* REFERENCE_PARTICLE_H_ */
