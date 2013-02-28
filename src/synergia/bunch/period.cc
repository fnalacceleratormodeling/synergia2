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

void apply_zcut(Bunch & bunch, double length)
{
           Bunch::State state=bunch.get_state(); 
           bunch.convert_to_state(Bunch::fixed_t_lab); 
           double half_length=0.5*length;	
    
            MArray2d_ref particles(bunch.get_local_particles());
            int discarded = 0;
            int local_num = bunch.get_local_num();
            for (int part = 0; part < local_num; ++part) {
                bool try_discard = true;
                while (try_discard) {
                    if  (fabs((particles[part][Bunch::z]))>half_length){
                          ++discarded;
                         --local_num;       
                         if (part == local_num) {
                        // No more particles left
                        try_discard = false;
                    } else {
                        // Move the last particle into this newly empty position
                        int last = local_num;
                        particles[part][0] = particles[last][0];
                        particles[part][1] = particles[last][1];
                        particles[part][2] = particles[last][2];
                        particles[part][3] = particles[last][3];
                        particles[part][4] = particles[last][4];
                        particles[part][5] = particles[last][5];
                        particles[part][6] = particles[last][6];
                    }       
                    } else{
                        try_discard = false;       
                    }
                                            
                } 
            }            
            bunch.set_local_num(local_num); 
            bunch.update_total_num();             
            bunch.convert_to_state(state); 
}
