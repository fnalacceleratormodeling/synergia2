#include <iostream>
#include "step.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/collective/impedance.h"
#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/get_pointer.hpp>
#include "synergia/collective/impedance.h"
#include "synergia/bunch/period.h"


Step::Step(double length) :
    operators(), time_fractions(), length(length)
{

}

void
Step::append(Operator_sptr operator_sptr, double time_fraction)
{
    operators.push_back(operator_sptr);
    time_fractions.push_back(time_fraction);
}

void
Step::append(Operators const& the_operators, double time_fraction)
{
    Operators tmp(the_operators);
    operators.splice(operators.end(), tmp);
    for (Operators::const_iterator it = the_operators.begin(); it
            != the_operators.end(); ++it) {
        time_fractions.push_back(time_fraction);
    }
}

void
Step::apply(Bunch & bunch)
{ 
    
   // int rank = Commxx().get_rank();
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
       //  if (rank==0) std::cout<<" operator name="<<(*it)->get_name()<<std::endl;
        // time [s] in accelerator frame
        double time = length / (bunch.get_reference_particle().get_beta()
                * pconstants::c);
        if ((*it)->get_name()=="impedance") {
            MArray1d bunch_means=Diagnostics::calculate_mean(bunch);
            Bunch_means bi;
            std::vector<Bunch_means> vbi;
            bi.x_mean=bunch_means[0];
            bi.y_mean=bunch_means[2];
            bi.z_mean=bunch_means[4];
            bi.realnum=bunch.get_real_num();           
            bi.bucket_index=bunch.get_bucket_index();
            vbi.push_back(bi);
            stored_vbunches.push_front(vbi);

            int nstored=(reinterpret_cast<Impedance*>(boost::get_pointer(*it)))->get_nstored_turns(); 
            if (stored_vbunches.size()>nstored) stored_vbunches.pop_back();
           
         } 
        (*it)->apply(bunch, (*fractions_it) * time, *this);
        if (bunch.is_z_periodic()){
            double plength=bunch.get_z_period_length();
            apply_longitudinal_periodicity(bunch, plength);
        }
        ++fractions_it;  
    }
}

void
Step::apply(Bunch & bunch, Multi_diagnostics & diagnostics)
{ 
    
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
         // std::cout<<" operator name="<<(*it)->get_name()<<std::endl;
        // time [s] in accelerator frame
        double time = length / (bunch.get_reference_particle().get_beta()
                * pconstants::c);
        if ((*it)->get_name()=="impedance") {
            MArray1d bunch_means=Diagnostics::calculate_mean(bunch);
            Bunch_means bi;
            std::vector<Bunch_means> vbi;
            bi.x_mean=bunch_means[0];
            bi.y_mean=bunch_means[2];
            bi.z_mean=bunch_means[4];
            bi.realnum=bunch.get_real_num();
            bi.bucket_index=bunch.get_bucket_index();
            vbi.push_back(bi);
            stored_vbunches.push_front(vbi);
           
            int nstored=(reinterpret_cast<Impedance*>(boost::get_pointer(*it)))->get_nstored_turns(); 
             if (stored_vbunches.size()>nstored) stored_vbunches.pop_back();
           
           // std::cout<<"name ="<< (*it)->get_name()<<" stored dim "<<stored_bunches.size()<<std::endl; 
           
         }
         for (Multi_diagnostics::iterator itd = diagnostics.begin(); itd
            != diagnostics.end(); ++itd) {
              
                 (*itd)->update_and_write();
          }
        (*it)->apply(bunch, (*fractions_it) * time, *this); 
         
         for (Multi_diagnostics::iterator itd = diagnostics.begin(); itd
            != diagnostics.end(); ++itd) {
              
                 (*itd)->update_and_write();
          }    
             
         if (bunch.is_z_periodic()){
            double plength=bunch.get_z_period_length();
            apply_longitudinal_periodicity(bunch, plength);
        }     
        ++fractions_it;
    }
}


void
Step::apply(Bunch_with_diagnostics_train & bunch_diag_train)
{ 
    int mrank = bunch_diag_train.get_master_comm().get_rank();
    int numbunches=bunch_diag_train.get_num_bunches();
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        if ((*it)->get_name()=="impedance") {
         Bunch_means bi;
         std::vector<Bunch_means> vbi_local(0);
         std::vector<Bunch_means> vbi(numbunches);
            for (int index = 0; index < numbunches; ++index) {
                if (bunch_diag_train.is_on_this_rank(index)) {
                    Bunch_sptr bunch_sptr=bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
                    MArray1d bunch_means=Diagnostics::calculate_mean(*bunch_sptr);
                    bi.x_mean=bunch_means[0];
                    bi.y_mean=bunch_means[2];
                    bi.z_mean=bunch_means[4];
                    bi.realnum=bunch_sptr->get_real_num();
                    bi.bucket_index=bunch_sptr->get_bucket_index();                
                    if  (bunch_diag_train.get_comm(index).get_rank() ==0)   vbi_local.push_back(bi);
                    ///only the rank 0 of every communicator sends the bi to all others
                }                
            }
           
            MPI_Datatype Bunch_means_type;
            MPI_Aint lb, extent;
            MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent); 
            MPI_Datatype type[2] = {MPI_DOUBLE, MPI_INT};
            int blocklen[2] = {4,1};
            MPI_Aint disp[2];
            disp[0]=0;
            disp[1]=4*extent;
            MPI_Type_create_struct(2,blocklen, disp, type, &Bunch_means_type);
            MPI_Type_commit(&Bunch_means_type); 
            
            Commxx master_comm=bunch_diag_train.get_master_comm();
          
            std::vector<int > counts=bunch_diag_train.get_proc_counts();
            std::vector<int > offsets=bunch_diag_train.get_proc_offsets();
            int error = MPI_Allgatherv(reinterpret_cast<void*>(&vbi_local[0]), vbi_local.size(), Bunch_means_type,  
                                         reinterpret_cast<void*>(&vbi[0]), &counts[0], &offsets[0], 
                                         Bunch_means_type, master_comm.get());
                       if (error != MPI_SUCCESS) {
                         throw std::runtime_error("step.cc::apply: MPI error in MPI_Allgatherv");
                       } 
                       
            MPI_Type_free(&Bunch_means_type);
            stored_vbunches.push_front(vbi);

            int nstored=(reinterpret_cast<Impedance*>(boost::get_pointer(*it)))->get_nstored_turns(); 
            if (stored_vbunches.size()>nstored) stored_vbunches.pop_back();
        

            
        }
                   
          
        
       for (int index = 0; index < numbunches; ++index) {
           if (bunch_diag_train.is_on_this_rank(index)) {
               Bunch_sptr bunch_sptr=bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
        // time [s] in accelerator frame
               double time = length / (bunch_sptr->get_reference_particle().get_beta() * pconstants::c);
               (*it)->apply(*bunch_sptr, (*fractions_it) * time, *this);
           }
       }
   
        for (int index = 0; index < numbunches; ++index) {
            if (bunch_diag_train.is_on_this_rank(index)) {   
                Bunch_sptr bunch_sptr=bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
                if (bunch_sptr->is_z_periodic()){
                    double plength=bunch_sptr->get_z_period_length();
                    apply_longitudinal_periodicity(*bunch_sptr, plength);
                } 
            }
         }    
      ++fractions_it;
  }  
}


Operators const&
Step::get_operators() const
{
    return operators;
}

std::list<double > const&
Step::get_time_fractions() const
{
    return time_fractions;
}

double
Step::get_length() const
{
    return length;
}

               
std::list< std::vector<Bunch_means> > const&
Step::get_stored_vbunches() const
{
  return  this->stored_vbunches;          
}



void
Step::print(int index) const
{
    std::cout << "step " << index << ":\n";
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        (*it)->print();
    }
}
