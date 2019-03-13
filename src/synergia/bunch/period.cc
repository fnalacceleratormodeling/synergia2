#include "period.h"

void apply_longitudinal_periodicity (Bunch & bunch, double length)
{
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
}

void apply_longitudinal_bucket_barrier(Bunch & bunch, double length)
{
}

void apply_zcut(Bunch & bunch, double length, Diagnostics_loss_sptr diag_loss_sptr)
{
           bool have_diagnostics=diag_loss_sptr.get();
                     
          if (have_diagnostics){
              if(bunch.is_bucket_index_assigned()){  
                  if (bunch.get_bucket_index() !=diag_loss_sptr->get_bunch().get_bucket_index())
                      throw std::runtime_error("apply_z_cut: diagnostics loss:  bunch bucket index difference"); 
               }
          }
             
           
           Bunch::State state=bunch.get_state(); 
           bunch.convert_to_state(Bunch::fixed_z_lab); 
           double beta=bunch. get_reference_particle().get_beta();
           double half_length=0.5*length/beta;	

           int b_index=-1; // AM: this value is written in the aperture_loss file when the bunch has no bucket index assigned
           if (bunch.is_bucket_index_assigned()) b_index=bunch.get_bucket_index();
           int repetition=bunch.get_reference_particle().get_repetition();
           double s=bunch.get_reference_particle().get_s();
           double s_n=bunch.get_reference_particle().get_s_n();
           MArray1d coords(boost::extents[7]);
           
    
            MArray2d_ref particles(bunch.get_local_particles());
            int discarded = 0;
            int local_num = bunch.get_local_num();
            for (int part = 0; part < local_num; ++part) {
                bool try_discard = true;
                while (try_discard) {
                    if  (std::abs((particles[part][Bunch::cdt]))>half_length){
                          ++discarded;
                         --local_num;  
                          if  (have_diagnostics){
                              coords[0]=particles[part][Bunch::x];
                              coords[1]=particles[part][Bunch::xp];
                              coords[2]=particles[part][Bunch::y];
                              coords[3]=particles[part][Bunch::yp];
                              coords[4]=particles[part][Bunch::z];
                              coords[5]=particles[part][Bunch::zp]; 
                              coords[6]=particles[part][Bunch::id];
                              diag_loss_sptr->update( b_index, repetition, s,s_n, coords ); 
                          }
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


