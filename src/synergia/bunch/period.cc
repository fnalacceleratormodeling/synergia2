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

    struct alg_z_period
    {
        Particles p;
        ConstParticleMasks masks;

        double length_cdt;
        double half_length;
        
        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                double tmp = p(i, 4) + half_length;
                if(tmp>0) p(i, 4) = fmod(tmp, length_cdt) - half_length;
                else      p(i, 4) = fmod(tmp, length_cdt) + half_length;
            }
        }
    };

}




void apply_longitudinal_periodicity(Bunch& bunch, double length)
{
    //Bunch::State state=bunch.get_state(); 
    //bunch.convert_to_state(Bunch::fixed_z_lab); 

    double beta = bunch.get_reference_particle().get_beta();      
    double length_cdt = length / beta;
    double half_length = 0.5 * length_cdt;

    if (bunch.get_local_num(ParticleGroup::regular))
    {
        alg_z_period alg{
            bunch.get_local_particles(ParticleGroup::regular),
            bunch.get_local_particle_masks(ParticleGroup::regular),
            length_cdt, half_length
        };

        Kokkos::parallel_for(
                bunch.get_local_num(ParticleGroup::regular), alg);
    }

    if (bunch.get_local_num(ParticleGroup::spectator))
    {
        alg_z_period alg{
            bunch.get_local_particles(ParticleGroup::spectator),
            bunch.get_local_particle_masks(ParticleGroup::spectator),
            length_cdt, half_length
        };

        Kokkos::parallel_for(
                bunch.get_local_num(ParticleGroup::spectator), alg);
    }

    //bunch.convert_to_state(state); 
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


