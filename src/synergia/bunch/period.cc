#include "period.h"
#include "synergia/bunch/bunch.h"

namespace
{
    struct Zcut_aperture
    {
        double half_length;

        KOKKOS_INLINE_FUNCTION
        bool discard(ConstParticles const& parts, int p) const
        {
           if  (std::abs(parts(p, 4)) > half_length) return true;
           else return false;
        }
    };
}




void apply_longitudinal_periodicity(Bunch& bunch, double length)
{
#if 0
       Bunch::State state=bunch.get_state(); 
       bunch.convert_to_state(Bunch::fixed_z_lab); 
       double beta=bunch. get_reference_particle().get_beta();      
       double  length_cdt=length/beta;
       double half_length=0.5*length_cdt;

        MArray2d_ref particles(bunch.get_local_particles());
        int local_num = bunch.get_local_num();
        for (int part = 0; part < local_num; ++part) {
            double tmp = particles[part][Bunch::cdt]+ half_length;
            if (tmp > 0) {
                particles[part][Bunch::cdt] = fmod(tmp, length_cdt) -half_length ;
            } else {
                particles[part][Bunch::cdt] = fmod(tmp, length_cdt) + half_length;
            }
        }
        bunch.convert_to_state(state); 
#endif
}

void apply_longitudinal_bucket_barrier(Bunch& bunch, double length)
{
}


void apply_zcut(Bunch& bunch, double length)
{
    //Bunch::State state = bunch.get_state(); 
    //bunch.convert_to_state(Bunch::fixed_z_lab); 
    double beta = bunch.get_reference_particle().get_beta();
    double half_length = 0.5 * length / beta;	
    bunch.apply_zcut(Zcut_aperture{half_length});
    //bunch.convert_to_state(state); 
}


