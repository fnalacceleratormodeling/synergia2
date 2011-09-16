#include "period.h"

void apply_longitudinal_periodicity (Bunch & bunch, double length)
{
       Bunch::State state=bunch.get_state(); 
       bunch.convert_to_state(Bunch::fixed_t_lab); 
       double half_length=0.5*length;	

        MArray2d_ref particles(bunch.get_local_particles());
        int local_num = bunch.get_local_num();
        for (int part = 0; part < local_num; ++part) {
            double tmp = particles[part][4]+ half_length;
            if (tmp > 0) {
                particles[part][4] = fmod(tmp, length) -half_length ;
            } else {
                particles[part][4] = fmod(tmp, length) + half_length;
            }
        }
        bunch.convert_to_state(state); 
}

