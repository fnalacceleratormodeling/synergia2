#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/multi_array_print.h"
#include <boost/archive/archive_exception.hpp>
#include "synergia/utils/multi_array_to_string.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/split_operator_stepper_choice.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/simulation/populate_stationary.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_phase_space_density.h"
#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/collective/impedance.h"
//#include "synergia/simulation/diagnostics_aperture.h"
#include "booster_options.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{ 
  
   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);     
  //  std::cout<<"run start on rank: "<<rank<<std::endl; 
    
    Options opts;
    opts.print();
  
    
    Lattice_sptr lattice_sptr(new Lattice());
    try {
        xml_load(*lattice_sptr, opts.lattice_file);
    }
    catch (std::exception& ex) {
        std::cerr <<"exception what "<<ex.what()<<std::endl;
        std::cerr << "cxx_booster: failed to find "<<opts.lattice_file<<"\n";
        std::cerr << "Run booster_xml.py to generate "<<opts.lattice_file<<"\n";
        return;
    }
    if (rank==0) { 
      std::cout<<" lattice file ="<<opts.lattice_file<<std::endl;
    }
  
   
    MPI_Barrier(MPI_COMM_WORLD);   
    double tini = MPI_Wtime();
    
    
    Lattice_elements quad_correctors_h;
    Lattice_elements quad_correctors_v;
    Lattice_elements sextupole_correctors_h;
    Lattice_elements sextupole_correctors_v;
   
    for (Lattice_elements::const_iterator it =
      lattice_sptr->get_elements().begin();
    it != lattice_sptr->get_elements().end(); ++it) {

          if (opts.chef_propagate) {
            (*it)->set_string_attribute("extractor_type", "chef_propagate");
          }
          else if (opts.chef_map) {
            (*it)->set_string_attribute("extractor_type", "chef_map");
          }
          
          std::string element_name=(*it)->get_name();
          std::string rf_name(element_name.begin(),element_name.begin()+2);
          if  (rf_name=="rf"){
            //std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
            (*it)->set_double_attribute("volt", opts.rf_voltage);
            (*it)->set_double_attribute("lag", 0.);
            (*it)->set_double_attribute("harmon", opts.harmon);
           // (*it)->print();
          }
          
          
          
          
          
          if (opts.if_aperture) {                          
            std::string mag_name(element_name.begin(),element_name.begin()+4);           
            if  (mag_name=="fmag") {               
              // std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
              (*it)->set_string_attribute("aperture_type","rectangular");
              (*it)->set_double_attribute("rectangular_aperture_width", 0.2);
              (*it)->set_double_attribute("rectangular_aperture_height", 2.*opts.aperture_f);
            }
            else if(mag_name=="dmag"){
              // std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
              (*it)->set_string_attribute("aperture_type","rectangular");
              (*it)->set_double_attribute("rectangular_aperture_width", 0.2);
              (*it)->set_double_attribute("rectangular_aperture_height", 2.*opts.aperture_d);
            }
            else{               
                (*it)->set_string_attribute("aperture_type","circular");
                (*it)->set_double_attribute("circular_aperture_radius", opts.aperture_l);              
            }            
          }//aperture
 
         std::string sex_name(element_name.begin(),element_name.begin()+3);         
          if ((sex_name=="sxs") &&  ((*it)->get_type()=="sextupole") ){            
          //  std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
           sextupole_correctors_h.push_back((*it));
          }
          if ( (sex_name=="sxl")   && ((*it)->get_type()=="sextupole") ){   
           //   std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
           sextupole_correctors_v.push_back((*it));
          }
          
          std::string quad_name(element_name.begin(),element_name.begin()+2); 
          std::string quad_namel(element_name.begin(),element_name.begin()+3); 
          if ( (quad_name=="ql")   && ((*it)->get_type()=="quadrupole") ) {
            // std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
            // (*it)->print();
             quad_correctors_h.push_back((*it));
          }      
            if ( (quad_name=="qs")   && ((*it)->get_type()=="quadrupole") && (quad_namel!="qss") && (quad_namel!="qsl") ) { 
            // std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
              (*it)->print();
             quad_correctors_v.push_back((*it));
          }
  
      
     }// for loop lattice elements

     if (rank==0) { 
          std::cout<<" number of horizontal tune quads correctors ="<<quad_correctors_h.size()<<std::endl;
          std::cout<<" number of vertical tune quads correctors ="<<quad_correctors_v.size()<<std::endl;   
          std::cout<<" number of horizontal chromaticity sextupole correctors ="<<sextupole_correctors_h.size()<<std::endl;
          std::cout<<" number of vertical chromaticity sextupole correctors ="<<sextupole_correctors_v.size()<<std::endl;   
     }
     
    
           
     Reference_particle reference_particle=lattice_sptr->get_reference_particle();        
     double lattice_length=lattice_sptr->get_length();       
     double beta = reference_particle.get_beta();
     double gamma = reference_particle.get_gamma();
     double energy = reference_particle.get_total_energy();
    
           
      if (rank==0) { 
         std::cout<<"    ***     before stepper   ***     "<<std::endl;    
         std::cout<<" beta="<<beta<<std::endl;
         std::cout<<" gamma="<<gamma<<std::endl;
         std::cout<<" lattice_length="<<lattice_length<<std::endl;
         std::cout<<" energy="<<energy<<std::endl;
         std::cout<<std::endl;
         std::cout<<"    ***    ***********    ***     "<<std::endl;
       }
      
    
      
    
    Stepper_sptr stepper_sptr;
    
    Dummy_collective_operator_sptr bpm_measure_sptr;
     if (opts.bpms){
        bpm_measure_sptr=Dummy_collective_operator_sptr(new Dummy_collective_operator("bmp_measure"));    
     }
    
    
    
     Space_charge_rectangular_sptr spc_f_sptr;
     Space_charge_rectangular_sptr spc_d_sptr;
     Space_charge_rectangular_sptr spc_l_sptr;     
     if (opts.space_charge){
      std::vector<int> grid_shape_f(opts.scgrid);
      double radiusx_f=0.1;
      double radiusy_f= opts.aperture_f;
      grid_shape_f[0] *= int(radiusx_f/opts.grid_ref_distance);
      grid_shape_f[1] *= int(radiusy_f/opts.grid_ref_distance);
      
      std::vector<double> pipe_size_f(3);
      pipe_size_f[0]=2.*radiusx_f;
      pipe_size_f[1]=2.*radiusy_f;
      pipe_size_f[2]=lattice_sptr->get_length()/opts.harmon;
      spc_f_sptr=Space_charge_rectangular_sptr(new Space_charge_rectangular(pipe_size_f, grid_shape_f));
      if (rank==0){
          std::cout<<std::endl;
          std::cout<<"pipe_size F magnet=["<<spc_f_sptr->get_pipe_size()[0]<<", "
                                           <<spc_f_sptr->get_pipe_size()[1]<<", "
                                           <<spc_f_sptr->get_pipe_size()[2]<<"]"
                                           <<std::endl;
          
         std::cout<<"grid for spc F magnet=["<< spc_f_sptr->get_grid_shape()[0]<<", "
                                           << spc_f_sptr->get_grid_shape()[1]<<", "
                                           << spc_f_sptr->get_grid_shape()[2]<<"]"         
                                            <<std::endl;
          std::cout<<"___________________________________________________________"<<std::endl;
      }  
      
      std::vector<int> grid_shape_d(opts.scgrid);
      double radiusx_d=0.1;
      double radiusy_d= opts.aperture_d;
      grid_shape_d[0] *= int(radiusx_d/opts.grid_ref_distance);
      grid_shape_d[1] *= int(radiusy_d/opts.grid_ref_distance);
      
      std::vector<double> pipe_size_d(3);
      pipe_size_d[0]=2.*radiusx_d;
      pipe_size_d[1]=2.*radiusy_d;
      //pipe_size_d[2]=lattice_simulator.get_bucket_length();
      pipe_size_d[2]=lattice_sptr->get_length()/opts.harmon;
      spc_d_sptr=Space_charge_rectangular_sptr(new Space_charge_rectangular(pipe_size_d, grid_shape_d));
      if (rank==0){
          std::cout<<std::endl;
          std::cout<<"pipe_size  D magnet=["<<spc_d_sptr->get_pipe_size()[0]<<", "
                                           <<spc_d_sptr->get_pipe_size()[1]<<", "
                                           <<spc_d_sptr->get_pipe_size()[2]<<"]"
                                           <<std::endl;
          
         std::cout<<"grid for spc D magnet=["<< spc_d_sptr->get_grid_shape()[0]<<", "
                                           << spc_d_sptr->get_grid_shape()[1]<<", "
                                           << spc_d_sptr->get_grid_shape()[2]<<"]"         
                                            <<std::endl;
          std::cout<<"___________________________________________________________"<<std::endl;
      }  

      std::vector<int> grid_shape_l(opts.scgrid_l);
      double radius_l=opts.aperture_l;
      std::vector<double> pipe_size_l(3);
      pipe_size_l[0]=2.*radius_l;
      pipe_size_l[1]=2.*radius_l;
      pipe_size_l[2]= lattice_sptr->get_length()/opts.harmon;
      spc_l_sptr=Space_charge_rectangular_sptr(new Space_charge_rectangular(pipe_size_l, grid_shape_l));
      if (rank==0){
          std::cout<<std::endl;
          std::cout<<"pipe_size L section=["<<spc_l_sptr->get_pipe_size()[0]<<", "
                                           <<spc_l_sptr->get_pipe_size()[1]<<", "
                                           <<spc_l_sptr->get_pipe_size()[2]<<"]"
                                           <<std::endl;
          
         std::cout<<"grid for spcL section=["<< spc_l_sptr->get_grid_shape()[0]<<", "
                                           << spc_l_sptr->get_grid_shape()[1]<<", "
                                           << spc_l_sptr->get_grid_shape()[2]<<"]"         
                                            <<std::endl;
          std::cout<<"___________________________________________________________"<<std::endl;
      }  
              
  }
  
 
  
  Impedance_sptr imped_f_sptr;
  Impedance_sptr imped_d_sptr;
  if (opts.impedance){
    int zgrid=1000;
    imped_f_sptr=Impedance_sptr(new Impedance(opts.wakefile_f, opts.waketype, zgrid, lattice_sptr->get_length(),
                                              lattice_sptr->get_length()/opts.harmon,opts.registred_turns,opts.full_machine));// opts.full_machine,wn));
                                              
    /*imped_f_sptr->get_wake_field_sptr()->multiply_xw_lead(0.);
    imped_f_sptr->get_wake_field_sptr()->multiply_xw_trail(0);
    imped_f_sptr->get_wake_field_sptr()->multiply_yw_lead(0);
    imped_f_sptr->get_wake_field_sptr()->multiply_yw_trail(0);
    imped_f_sptr->get_wake_field_sptr()->multiply_z_wake(0);      */                                    
                                             
     if (rank==0){
         std::cout<<std::endl;
         std::cout<<"WAKES FOR F MAGNET read from "<<imped_f_sptr->get_wake_field_sptr()->get_wake_file_name()<<std::endl;
         std::cout<<"F mag orbith length="<<imped_f_sptr->get_orbit_length()<<std::endl;
         std::cout<<"F mag z_grid="<<imped_f_sptr->get_z_grid()<<std::endl;
         std::cout<<"F mag stored turns="<<imped_f_sptr->get_nstored_turns()<<std::endl;
     }
  
     imped_d_sptr=Impedance_sptr(new Impedance(opts.wakefile_d, opts.waketype, zgrid, lattice_sptr->get_length(),
                                               lattice_sptr->get_length()/opts.harmon,opts.registred_turns,opts.full_machine));// opts.full_machine,wn));
  
     /*                                           
    imped_d_sptr->get_wake_field_sptr()->multiply_xw_lead(0.); 
    imped_d_sptr->get_wake_field_sptr()->multiply_xw_trail(0);
    imped_d_sptr->get_wake_field_sptr()->multiply_yw_lead(0);
    imped_d_sptr->get_wake_field_sptr()->multiply_yw_trail(0);
    imped_d_sptr->get_wake_field_sptr()->multiply_z_wake(0);  */  
     
     if (rank==0){
        std::cout<<std::endl;
        std::cout<<"WAKES FOR D MAGNET read from "<<imped_d_sptr->get_wake_field_sptr()->get_wake_file_name()<<std::endl;
        std::cout<<"D mag orbith length="<<imped_d_sptr->get_orbit_length()<<std::endl;
        std::cout<<"D mag z_grid="<<imped_d_sptr->get_z_grid()<<std::endl;
        std::cout<<"D mag stored turns="<<imped_d_sptr->get_nstored_turns()<<std::endl;
    }                                          
                                              
  }
   
//    Dummy_collective_operator_sptr no_op_sptr;
   if (  (opts.space_charge) ||  (opts.impedance) || (opts.bpms)  ) {
      List_choice_map list_choice_map;
 
      Collective_operators operators_fmag;
      if (opts.impedance) operators_fmag.push_back(imped_f_sptr);
      if (opts.space_charge) operators_fmag.push_back(spc_f_sptr);
      Kicks  kicks_fmag(operators_fmag, opts.steps_per_fmag);
      
      Collective_operators operators_dmag;
      if (opts.impedance) operators_dmag.push_back(imped_d_sptr);
      if (opts.space_charge) operators_dmag.push_back(spc_d_sptr);
      Kicks  kicks_dmag(operators_dmag, opts.steps_per_dmag);
      

      int steps_per_bpm=1;
      Collective_operators operators_bpm;
      if (opts.bpms) operators_bpm.push_back(bpm_measure_sptr);
      Kicks  kicks_bpm(operators_bpm, steps_per_bpm);
      
      
      for (Lattice_elements::const_iterator latt_it =
          lattice_sptr->get_elements().begin();
          latt_it !=lattice_sptr->get_elements().end(); ++latt_it){
            std::string element_name=(*latt_it)->get_name();  
            if ( (opts.space_charge) || (opts.impedance) ) {             
              if (element_name.substr(0,4)=="fmag") list_choice_map[element_name]=kicks_fmag;                 
              if (element_name.substr(0,4)=="dmag") list_choice_map[element_name]=kicks_dmag;                  
            }
            if (opts.bpms) {
              if(element_name.substr(0,3)=="bpm") {
                  list_choice_map[element_name]=kicks_bpm;
                //  if (rank==0) std::cout<<" BPMS are: "<<(*latt_it)->get_name()<<std::endl;
              }
            }
      }
      
      Collective_operators operators_else;
      if (opts.space_charge) operators_else.push_back(spc_l_sptr);
      Kicks  kicks_else(operators_else, opts.num_steps_else);
      list_choice_map["else"]=kicks_else;
      
      if (rank==0){
          std::cout<<std::endl; 
          std::cout<<"Split_operator_stepper_choice created"<<std::endl;
          if (opts.bpms) std::cout<<"measurements at BPMS elements"<<std::endl;
          if ( (opts.space_charge) || (opts.impedance) ) {
            std::cout<<"steps_per_fmag="<<opts.steps_per_fmag<<std::endl;
            std::cout<<"steps_per_dmag="<<opts.steps_per_dmag<<std::endl;
            std::cout<<"num_steps_else="<<opts.num_steps_else<<std::endl; 
          }
          std::cout<<"___________________________________________________________"<<std::endl;
      }         
      stepper_sptr=Stepper_sptr(new Split_operator_stepper_choice(lattice_sptr, opts.map_order, list_choice_map, true));
  }
  else{ 
   // no_op_sptr=Dummy_collective_operator_sptr(new Dummy_collective_operator("stub"));
   // int numst=1;    
   // stepper_sptr=Stepper_sptr(new Split_operator_stepper_elements(lattice_sptr, opts.map_order,  no_op_sptr,  numst ));   
     stepper_sptr=Stepper_sptr(new Independent_stepper(lattice_sptr, opts.map_order, opts.num_steps));  
  
                    
      if (rank==0) {
          std::cout<<"no collective effects, no bpm, stepper element operator"<<std::endl;
          std::cout<<"no collective effects, no bpm, independent operator"<<std::endl;
          std::cout<<"number of steps="<<stepper_sptr->get_steps().size()<<std::endl;
        }
    }
    
   //  stepper_sptr->get_lattice_simulator().print_lattice_functions();
     stepper_sptr->get_lattice_simulator().register_closed_orbit();             
     stepper_sptr->get_lattice_simulator().set_rf_bucket_length();
  
     if (rank==0) { 
       std::cout<<"stepper:lattice_simulator: map order="<< stepper_sptr->get_lattice_simulator().get_map_order() <<std::endl;
    // stepper_sptr->print();
       std::cout<<"stepper rf frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()<<std::endl;
     }
  
 
     if (opts.adjust_tunes){
          if (rank==0) std::cout<<"adjusting tunes"<<std::endl;
          stepper_sptr->get_lattice_simulator().adjust_tunes_chef(opts.tune_h, opts.tune_v,
                                                quad_correctors_h,  quad_correctors_v, 10, 1e-4, 3);  
          if (rank==0)  std::cout<<"tunes adjusted"<<std::endl; 
     }
 
     if (opts.adjust_chromaticity){
         if (rank==0) std::cout<<"adjusting chromaticity"<<std::endl;
         stepper_sptr->get_lattice_simulator().adjust_chromaticities(opts.chrom_h, opts.chrom_v, 
                                                     sextupole_correctors_h,sextupole_correctors_v);
         if (rank==0)  std::cout<<"chromaticity adjusted"<<std::endl;                                           
    }
 
   
          
     stepper_sptr->print_cs_step_betas();
     reference_particle=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_reference_particle();
     lattice_length=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_length();       
     beta = reference_particle.get_beta();
     gamma = reference_particle.get_gamma();
     energy = reference_particle.get_total_energy();
     double bunch_sp=stepper_sptr->get_lattice_simulator().get_bucket_length();
     MArray1d clo=stepper_sptr->get_lattice_simulator().get_closed_orbit();
     std::vector<double > actions= stepper_sptr->get_lattice_simulator().get_stationary_actions(opts.xrms, opts.yrms, opts.zrms/beta);
      if ((actions[0]<=0.) || (actions[1]<=0.) || (actions[2]<=0.)) throw
               std::runtime_error("get_stationary_actions can't satisfy requested moments");                    
      if (rank==0) { 
        std::cout<<"    ********    stepper  lattice  ************     "<<std::endl;
        std::cout<<std::endl;
        multi_array_print(clo,"closed_orbit");  
        std::cout<<" bunch spacing="<<bunch_sp<<std::endl;
        std::cout<<" rf bucket length="<<stepper_sptr->get_lattice_simulator().get_bucket_length()<<std::endl;        
        std::cout<<" rf frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()/1.e6<<" MHz"<<std::endl;     
        std::cout<<" beta="<<beta<<std::endl;
        std::cout<<" gamma="<<gamma<<std::endl;
        std::cout<<" lattice_length="<<lattice_length<<std::endl;
        std::cout<<" closed_orbit_length="<<stepper_sptr->get_lattice_simulator().get_closed_orbit_length()<<std::endl;
        std::cout<<" energy="<<energy<<std::endl;
        std::cout<<" actions= ("<< actions[0]<<", "<< actions[1]<<", "<<actions[2]<<")"<<std::endl;
        std::cout<<std::endl;
        std::cout<<"    ***********************************     "<<std::endl;
        std::cout<<std::endl;
       }



      MArray2d one_turn_map=stepper_sptr->get_lattice_simulator().get_linear_one_turn_map();
      MArray2d correlation_matrix=get_correlation_matrix(one_turn_map,opts.xrms, opts.yrms, opts.zrms, beta);
      if (rank==0){
         std::cout<<" correlation matrix="<<multi_array_to_string(correlation_matrix)<<std::endl;
       }
      int map_size=one_turn_map.size();
      Eigen::MatrixXd eigen_map(map_size,map_size);
      for (int i=0;i<map_size;++i){
        for (int j=0;j<map_size;++j){
          eigen_map(i,j)=one_turn_map[i][j];   
        }
      }
      if (rank==0){             
            std::cout<<" one turn map="<<multi_array_to_string(one_turn_map)<<std::endl;
            std::cout<<"eigenvalues of the one_turn_map:"
            <<"\n"<<"**********************"<<std::endl;
            for (int i=0;i<map_size;++i){
              std::complex<double> vv=eigen_map.eigenvalues()(i);
              std::cout<<"eigenvalue["<<i<<"]="<<vv
              <<",  absolute value="<<abs(vv)
              <<",  fractional tune="<<fabs(log(vv).imag()/(2.*mconstants::pi))
              <<",  "<<1.-fabs(log(vv).imag()/(2.*mconstants::pi))<<std::endl;
            }
            std::cout<<"**********************"<<std::endl;
      }

    
      if (opts.tunes_and_chroms){        
          double chef_frac_tunex=stepper_sptr->get_lattice_simulator().get_horizontal_tune();
          double chef_frac_tuney=stepper_sptr->get_lattice_simulator().get_vertical_tune();      
          double chef_eigen_tunex=stepper_sptr->get_lattice_simulator().get_horizontal_tune(1);
          double chef_eigen_tuney=stepper_sptr->get_lattice_simulator().get_vertical_tune(1);        
          double horizontal_chromaticity=stepper_sptr->get_lattice_simulator().get_horizontal_chromaticity();
          double vertical_chromaticity=stepper_sptr->get_lattice_simulator().get_vertical_chromaticity();                
          double momentum_compaction=stepper_sptr->get_lattice_simulator().get_momentum_compaction();
          double slip_factor=stepper_sptr->get_lattice_simulator().get_slip_factor();  
      
      
          if (rank==0){              
              std::cout<< "chef FracTune x: "<< chef_frac_tunex<< ", EigenTune x: "<< chef_eigen_tunex<<std::endl; 
              std::cout<< "chef FracTune y: "<< chef_frac_tuney<< ", EigenTune y: "<< chef_eigen_tuney<<std::endl;                   
              std::cout<< "horizontal chromaticity: "<< horizontal_chromaticity<<std::endl;
              std::cout<< "vertical   chromaticity: "<< vertical_chromaticity<<std::endl;
              std::cout<< "momentum compaction: "<< momentum_compaction<<std::endl;
              std::cout<< "slip factor: "<< slip_factor<<std::endl;
              
              
            }   
 } 
  


  Bunches bunches; 
  Commxx_sptr parent_comm_sptr(new Commxx); 
  Commxxs comms(generate_subcomms(parent_comm_sptr, opts.num_bunches));
  for (int i=0;i<opts.num_bunches;++i){
     Commxx_sptr commx=comms[i];      
       Bunch_sptr bunch_sptr=Bunch_sptr(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
          opts.num_macroparticles, opts.num_real_particles, commx));  
       bunch_sptr->set_bucket_index(i);
       if (opts.bunch_periodic){           
          bunch_sptr->set_z_period_length(stepper_sptr->get_lattice_simulator().get_bucket_length()) ;                    
      }
      else{
          bunch_sptr->set_longitudinal_aperture_length(stepper_sptr->get_lattice_simulator().get_bucket_length()) ;
      }       
      if (commx->has_this_rank()){
          Random_distribution dist(opts.seed,*commx);
          MArray1d input_means(boost::extents[6]);
          for(int imean=0;imean<6;++imean){
                 input_means[imean]=0.;
          }
          populate_6d(dist, *bunch_sptr, input_means, correlation_matrix);
       
      }      
      bunches.push_back(bunch_sptr);             
  }   
 
  Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, stepper_sptr->get_lattice_simulator().get_bucket_length()));
  Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);
  
 
  
  std::list<int> step_numbers;
  int step_number=0;
    for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it != stepper_sptr->get_steps().end(); ++it) {
        ++step_number;      
        if (opts.bpms){
           for (Operators::const_iterator o_it = (*it)->get_operators().begin(); o_it
                != (*it)->get_operators().end(); ++o_it) {
                  if ((*o_it)->get_name()=="bmp_measure") step_numbers.push_back(step_number);
           }
        }
        else{
          step_numbers.push_back(step_number);
        }
    }
   if (rank==0){
     std::cout<<" number of measurements per turn="<<step_numbers.size()<<std::endl;
   }
   
// #if 0  
  for (int i=0;i<bunch_train_sptr->get_size();++i){
    std::stringstream bunch_label;
    bunch_label<<i;
    if (opts.save_bunch){
       bunch_label<<"_m"<<opts.map_order;
       Diagnostics_sptr initial_particles_sptr(new Diagnostics_particles("initial_bunch"+bunch_label.str()+".h5"));
       initial_particles_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
       initial_particles_sptr->update_and_write();     
    }
   
     Diagnostics_sptr diagnostics_full2_sptr(new Diagnostics_full2("step_full2_"+bunch_label.str()+".h5"));
     bunch_train_simulator.add_per_step(i,diagnostics_full2_sptr, step_numbers ); 
     
   
     Diagnostics_sptr diagnostics_particles_sptr(new Diagnostics_particles("turn_particles_"+bunch_label.str()+".h5"));
     bunch_train_simulator.add_per_turn(i,diagnostics_particles_sptr,opts.turn_period);
         
     if(opts.turn_track){
        Diagnostics_sptr diagnostics_bulk_track_sptr(new Diagnostics_bulk_track("bulk_track_"+bunch_label.str()+".h5",
                                                                                opts.num_macroparticles));
        bunch_train_simulator.add_per_turn(i,diagnostics_bulk_track_sptr);  
     }
     if(opts.phase_space){
        int grid_z=20;
        int grid_zp=20;
        double z_nsigma=4.0;
        double zp_nsigma=4.0;
        Diagnostics_sptr diagnostics_phase_space_density_sptr(new 
                  Diagnostics_phase_space_density("phase_space_density_"+bunch_label.str()+".h5",
                                                  grid_z,grid_zp,z_nsigma,zp_nsigma));                                               
        bunch_train_simulator.add_per_turn(i, diagnostics_phase_space_density_sptr);       
     }
    
     if(opts.space_charge){
        if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().has_this_rank()){
           Commxx_sptr comm_spc_sptr=
                   make_optimal_spc_comm(bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm_sptr(),
                           opts.spc_comm_size, opts.equally_spread );
            spc_f_sptr->set_fftw_helper(comm_spc_sptr,opts.equally_spread);
            spc_d_sptr->set_fftw_helper(comm_spc_sptr,opts.equally_spread);
            spc_l_sptr->set_fftw_helper(comm_spc_sptr,opts.equally_spread); 
             if (opts.spc_tuneshift){
               Diagnostics_space_charge_rectangular_sptr sp_diagnostics_sptr
                  (new Diagnostics_space_charge_rectangular("space_charge_diagnostics_"+bunch_label.str()+".h5"));
               sp_diagnostics_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
               spc_f_sptr->add_diagnostics(sp_diagnostics_sptr);
               spc_d_sptr->add_diagnostics(sp_diagnostics_sptr);
               spc_l_sptr->add_diagnostics(sp_diagnostics_sptr);
           }// opts.spc_tuneshift                         
        }// has rank       
     }// space charge

     //adjust means   
     if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().has_this_rank()){
        MArray1d bunch_means=Core_diagnostics::calculate_mean(*bunch_train_simulator.get_bunch_train().get_bunches()[i]); 
        MArray2d_ref particles( bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_local_particles());
        for (int part=0;part<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_local_num();++part){
            particles[part][0] +=-bunch_means[0]+opts.xoffset;//+clo[0];
            particles[part][1] +=-bunch_means[1];//+clo[1];
            particles[part][2] +=-bunch_means[2]+opts.yoffset;//+clo[2];
            particles[part][3] +=-bunch_means[3];//+clo[3];
            particles[part][4] +=-bunch_means[4]+opts.zoffset;//+beta*clo[4];
            particles[part][5] +=-bunch_means[5];//+clo[5];
        }
        bunch_means=Core_diagnostics::calculate_mean(*bunch_train_simulator.get_bunch_train().get_bunches()[i]); 
        MArray1d bunch_stds=Core_diagnostics::calculate_std(*bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means);
        MArray2d bunch_mom2=Core_diagnostics::calculate_mom2(*bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means);
        if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().get_rank()==0) {
            std::cout<<std::endl; 
            std::cout<<"bunch # "<<i<<" is perriodic="<<bunch_train_simulator.get_bunch_train().get_bunches()[i]->is_z_periodic()<<std::endl;
            std::cout<<"bunch # "<<i<<" has longitudinal aperture="<<bunch_train_simulator.get_bunch_train().get_bunches()[i]->has_longitudinal_aperture()<<std::endl;
            std::cout<<"bunch # "<<i<<" number of real  particles= "
                  <<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_real_num()<<std::endl;
            std::cout<<"bunch # "<<i<<" number of macroparticles= "
                  <<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_total_num()<<std::endl;
            std::cout<<"bunch # "<<i<<" bucket index= "
                  <<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_bucket_index()<<std::endl;
            std::cout<<"bunch # "<<i<<" initial offsets (x,xp,y,yp,ct,dpp)=("
                  <<bunch_means[0]<<", "<<bunch_means[1]<<", "<<bunch_means[2]<<", "
                  <<bunch_means[3]<<", "<<bunch_means[4]<<", "<<bunch_means[5]
                  <<") [meters]"<<std::endl;
              std::cout<<"bunch # "<<i<<" initial stds (xrms,yrms,zrms)=("
                  <<bunch_stds[0]<<", "<<bunch_stds[2]<<", "<<beta*bunch_stds[4]<<") [meters]"<<std::endl;                   
             std::cout<<"bunch # "<<i<<" :"<<std::endl;
             Core_diagnostics::print_bunch_parameters(bunch_mom2, beta);
             std::cout<<"___________________________________________________________"<<std::endl; 
            
        }
     } // has rank

  }// for i
  if ((opts.num_bunches>1) && (rank==0)) std::cout<<"train bunch space="<<bunch_train_simulator.get_bunch_train().get_spacings()[0]<<std::endl;
 
  Propagator propagator(stepper_sptr);
  propagator.set_checkpoint_period(opts.checkpointperiod);
  propagator.set_concurrent_io(opts.concurrentio);
  propagator.propagate(bunch_train_simulator,opts.num_turns, opts.maxturns, opts.verbosity);
//#endif 
    MPI_Barrier(MPI_COMM_WORLD); 
    double   tfinal = MPI_Wtime();
          if (rank==0){ 
             std::ios_base::fmtflags old_flags(std::cout.flags());
              std::cout << std::scientific << std::setprecision(8);
              std::cout<<" time total ="<<tfinal-tini<<std::endl;
              std::cout.flags(old_flags);
          }
}

int
main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   run();
   MPI_Finalize();
   return 0;
      
}
