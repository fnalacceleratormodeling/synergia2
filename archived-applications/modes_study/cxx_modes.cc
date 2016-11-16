#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/multi_array_to_string.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/split_operator_stepper_choice.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/populate_stationary.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_phase_space_density.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/collective/space_charge_2d_bassetti_erskine.h"
#include "synergia/collective/space_charge_rectangular.h"
//#include "synergia/collective/dummy_collective_alex.h"
//#include "synergia/simulation/dummy_collective_test.h"
#include "modes_options.h"
#include "excite_mode_actions.h"
//#include "diagnostics_deltae_particles.h"

void
run()
{

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  

  std::cout<<"run start on rank: "<<rank<<std::endl; 
 
    
  Options opts;
  opts.print();
  
  Lattice_sptr lattice_sptr(new Lattice());
  
 
  
  try {
      xml_load(*lattice_sptr, opts.lattice_file);
  }
  catch (std::exception& ex) {
        std::cerr <<" boost xml load error "<<ex.what()<<std::endl;
        std::cerr << "cxx_booster: failed to find "<<opts.lattice_file<<"\n";
        std::cerr << "Run modes_xml.py to generate modes_lattice.xml\n";
        return;
  }

  
  if (rank==0) std::cout<<" latice loaded"<<std::endl; 

  Reference_particle reference_particle=lattice_sptr->get_reference_particle()  ;  
  double lattice_length=lattice_sptr->get_length();     
  
  double beta = reference_particle.get_beta();
  double gamma = reference_particle.get_gamma();
  double energy = reference_particle.get_total_energy();
  double freq=opts.harmon*beta*pconstants::c/lattice_length;
  if (rank==0){
    std::cout<< "lattice length [m]: "<< lattice_length<<std::endl; 
    std::cout<< "reference particle total energy: "<< energy<<std::endl; 
    std::cout<< "reference particle beta: "<< beta<<std::endl; 
    std::cout<< "reference particle gamma: "<<gamma<<std::endl; 
  }
  
 
 
  for (Lattice_elements::const_iterator it =
                lattice_sptr->get_elements().begin();
                it != lattice_sptr->get_elements().end(); ++it) {
    
      if (opts.chef_propagate) {
        (*it)->set_string_attribute("extractor_type", "chef_propagate");
      }
      else if (opts.chef_map) {
        (*it)->set_string_attribute("extractor_type", "chef_map");
      }

      if ( (*it)->get_type()=="rfcavity") {
        (*it)->set_double_attribute("volt", opts.rf_voltage);
       // (*it)->set_double_attribute("freq", freq);
        (*it)->set_double_attribute("lag", 0.);
        (*it)->set_double_attribute("harmon", opts.harmon);
      }  
       
      if (opts.if_aperture) {        
         (*it)->set_string_attribute("aperture_type","circular");
         (*it)->set_double_attribute("circular_aperture_radius", opts.aperture);
      }

     if ( (*it)->get_type()=="octupole") {
          (*it)->set_string_attribute("extractor_type", "chef_propagate");
          if ( (*it)->get_name()=="e1") (*it)->set_double_attribute("k3",opts.k3_e1);
          //if ( (*it)->get_name()=="e2") (*it)->set_double_attribute("k3",opts.k3_e2);
          double k3_e2=-0.59968587*opts.k3_e1;      // for \Delta_Q=a(J_x+J_y). i.e. a=b
          if ( (*it)->get_name()=="e2") (*it)->set_double_attribute("k3",k3_e2);
          
          //std::cout<<(*it)->get_name()<<" : "<<(*it)->get_string_attribute("extractor_type")
         // <<"  k3="<<(*it)->get_double_attribute("k3") <<std::endl;
      }
      
  }
  
//    for (Lattice_elements::const_iterator it =
//                 lattice_sptr->get_elements().begin();
//                 it != lattice_sptr->get_elements().end(); ++it) {
//       if ( (*it)->get_type()=="octupole") {
//       std::cout<<(*it)->get_name()<<" : "<<(*it)->get_string_attribute("extractor_type")
//           <<"  k3begin="<<(*it)->get_double_attribute("k3") <<std::endl;
//       }
//    }
 
  Lattice_simulator lattice_simulator(lattice_sptr, opts.map_order);
  //lattice_simulator.print_lattice_functions();
  lattice_simulator.register_closed_orbit();
  


   Stepper_sptr stepper_sptr;
   
   Space_charge_2d_bassetti_erskine_sptr spcbe_sptr;
   Dummy_collective_operator_sptr spc2dh_sptr;
   //Space_charge_2d_open_hockney_sptr spc2dh_sptr;
   Space_charge_3d_open_hockney_sptr spc3dh_sptr;
   Space_charge_rectangular_sptr spcrec_sptr;
   
   if (opts.space_charge_be){
      spcbe_sptr=Space_charge_2d_bassetti_erskine_sptr( new Space_charge_2d_bassetti_erskine() );
      if (rank==0){
         std::cout<<std::endl;
         std::cout<< "BASSETTI-ERSKINE SPACE CHARGE SOLVER"<<std::endl;
         std::cout<<"___________________________________________________________"<<std::endl;    
      }
   }   
   else if (opts.space_charge_2dh){
     spc2dh_sptr = Dummy_collective_operator_sptr(new Dummy_collective_operator("dummilica"));
//       std::vector<int> grid_shape(opts.scgrid);
//       Commxx_sptr spc_comm_sptr(new Commxx);
//       spc2dh_sptr = Space_charge_2d_open_hockney_sptr( new Space_charge_2d_open_hockney(spc_comm_sptr, grid_shape ));
//       if (rank==0){
//          std::cout<<std::endl;
//          std::cout<< "2D HOCKNEY SOLVER"<<std::endl;
//          std::cout<<"grid for spc section=["<< grid_shape[0]<<", "
//                                             << grid_shape[1]<<", "
//                                             << grid_shape[2]<<"]"         
//                                               <<std::endl;
//          std::cout<<"___________________________________________________________"<<std::endl;                                   
//       }
   }   
    else if (opts.space_charge_3dh){
      std::vector<int> grid_shape(opts.scgrid);
      bool longitudinal_kicks = false;
      bool periodic_z = false;
      double z_period=gamma*lattice_length/opts.harmon;
      bool grid_entire_period = false;
      double nsigma=8.;
     //  spc3dh_sptr = Space_charge_3d_open_hockney_sptr( new Space_charge_3d_open_hockney(grid_shape));
       spc3dh_sptr = Space_charge_3d_open_hockney_sptr( 
       new Space_charge_3d_open_hockney(grid_shape,longitudinal_kicks,periodic_z,z_period,grid_entire_period,nsigma));
       spc3dh_sptr->set_green_fn_type(Space_charge_3d_open_hockney::linear);
       //spc3dh_sptr->set_green_fn_type(Space_charge_3d_open_hockney::pointlike);
//       spc3dh_sptr->set_kick_frame(Bunch::fixed_z_lab);
       // spc3dh_sptr->set_green_fn_type(Space_charge_3d_open_hockney::pointlike);
       //spc3dh_sptr->set_kick_frame(Bunch::fixed_t_bunch);
   
     
       if (rank==0){
         std::cout<<std::endl;
         std::cout<< "3D HOCKNEY SOLVER"<<std::endl;
         std::cout<<"grid for spc section=["<< grid_shape[0]<<", "
                                            << grid_shape[1]<<", "
                                            << grid_shape[2]<<"]"         
                                              <<std::endl;
         std::cout<<"___________________________________________________________"<<std::endl;                                   
      }
    }   
    else if (opts.space_charge_rec){
       
        std::vector<int> grid_shape(opts.scgrid);
        double radius=opts.aperture;
        std::vector<double> pipe_size(3);
        pipe_size[0]=2.*radius;
        pipe_size[1]=2.*radius;
        pipe_size[2]=lattice_simulator.get_bucket_length();
        
      
        
        
        spcrec_sptr=Space_charge_rectangular_sptr(new Space_charge_rectangular(pipe_size, grid_shape));      
        if (rank==0){
            std::cout<<std::endl;
            std::cout<< "RECTANGULAR SOLVER"<<std::endl;
            std::cout<<"pipe_size  section=["<<spcrec_sptr->get_pipe_size()[0]<<", "
                                            <<spcrec_sptr->get_pipe_size()[1]<<", "
                                            <<spcrec_sptr->get_pipe_size()[2]<<"]"
                                            <<std::endl;
            
            std::cout<<"grid for spc section=["<< spcrec_sptr->get_grid_shape()[0]<<", "
                                            << spcrec_sptr->get_grid_shape()[1]<<", "
                                            << spcrec_sptr->get_grid_shape()[2]<<"]"         
                                              <<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
      }          
  }

 
    Dummy_collective_operator_sptr dummy_sptr;   
    dummy_sptr=Dummy_collective_operator_sptr(new Dummy_collective_operator("dummy"));    
    

  bool octupole_dummy=true;
  if ((opts.space_charge_be) || (opts.space_charge_2dh) || (opts.space_charge_3dh) || (opts.space_charge_rec) || octupole_dummy){
      List_choice_map list_choice_map;
 
      Collective_operators operators_quad;
      if (opts.space_charge_be) operators_quad.push_back(spcbe_sptr);
      if (opts.space_charge_2dh) operators_quad.push_back(spc2dh_sptr);
      if (opts.space_charge_3dh) operators_quad.push_back(spc3dh_sptr);
      if (opts.space_charge_rec) operators_quad.push_back(spcrec_sptr);
      Kicks  kicks_quad(operators_quad,opts.steps_per_quad);
      
      Collective_operators operators_sbend;
      if (opts.space_charge_be) operators_sbend.push_back(spcbe_sptr);
      if (opts.space_charge_2dh) operators_quad.push_back(spc2dh_sptr);
      if (opts.space_charge_3dh) operators_quad.push_back(spc3dh_sptr);
      if (opts.space_charge_rec) operators_quad.push_back(spcrec_sptr);
      Kicks  kicks_sbend(operators_sbend,opts.steps_per_sbend);
      
      int steps_per_oct=1;
      Collective_operators operators_oct;
      operators_oct.push_back(dummy_sptr);
      Kicks  kicks_oct(operators_oct, steps_per_oct);
      
      
      for (Lattice_elements::const_iterator latt_it =
          lattice_sptr->get_elements().begin();
          latt_it !=lattice_sptr->get_elements().end(); ++latt_it){
        
          if (((*latt_it)->get_name()=="f") || ((*latt_it)->get_name()=="d")) {             
                 // std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
              list_choice_map[(*latt_it)->get_name()]=kicks_quad;
          }
          
          if ((*latt_it)->get_name()=="b") list_choice_map[(*latt_it)->get_name()]=kicks_sbend;
          
          if ( ((*latt_it)->get_name()=="e1") || ((*latt_it)->get_name()=="e2") ){
                 // std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
              list_choice_map[(*latt_it)->get_name()]=kicks_oct;
          }
      }
      
      Collective_operators operators_else;
      if (opts.space_charge_be) operators_else.push_back(spcbe_sptr);
      if (opts.space_charge_2dh) operators_else.push_back(spc2dh_sptr);
      if (opts.space_charge_3dh) operators_else.push_back(spc3dh_sptr);
      if (opts.space_charge_rec) operators_else.push_back(spcrec_sptr);
      Kicks  kicks_else(operators_else,opts.num_steps_else);
      list_choice_map["else"]=kicks_else;
      
      if (rank==0){
          std::cout<<std::endl; 
          std::cout<<"steps_per_quad="<<opts.steps_per_quad<<std::endl;
          std::cout<<"steps_per_sbend="<<opts.steps_per_sbend<<std::endl;
          std::cout<<"num_steps_else="<<opts.num_steps_else<<std::endl;
          std::cout<<"___________________________________________________________"<<std::endl;
      }
       
       
     stepper_sptr=Stepper_sptr(new Split_operator_stepper_choice(lattice_simulator,list_choice_map, true));      
    //  stepper_sptr= Stepper_sptr(new Split_operator_stepper(lattice_simulator,  spc3dh_sptr, 40));
   
   
   
  }
  else{
    stepper_sptr=Stepper_sptr(new Independent_stepper(lattice_simulator, opts.num_steps));
    if (rank==0) {
      std::cout<<"no collective effects, no bpm, independent operator"<<std::endl;
      std::cout<<"number of steps="<<stepper_sptr->get_steps().size()<<std::endl;
      std::cout<<"stepper:lattice_simulator: map order="<< stepper_sptr->get_lattice_simulator().get_map_order() <<std::endl;
    }
  }  

  //stepper_sptr->print();
  stepper_sptr->get_lattice_simulator().register_closed_orbit();
  double bunch_sp= stepper_sptr->get_lattice_simulator().get_bucket_length();
  if (rank==0){
      std::cout<< "RF frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()<<std::endl;
      std::cout<< "bucket length: "<< bunch_sp<<std::endl;
  }

 
  
    MArray2d one_turn_map=stepper_sptr->get_lattice_simulator().get_linear_one_turn_map();
    int map_size=one_turn_map.size();
    Eigen::MatrixXd eigen_map(map_size,map_size);
    for (int i=0;i<map_size;++i){
      for (int j=0;j<map_size;++j){
        eigen_map(i,j)=one_turn_map[i][j];   
      }
    }

    double  cosmu=0.5*(one_turn_map[0][0]+one_turn_map[1][1]); 
    if (fabs(cosmu)>1.) throw std::runtime_error(" alpha_x: cosmu larger than zero");
    double sinmu=sqrt(1.-cosmu*cosmu);
    double alpha_x=(one_turn_map[0][0]-cosmu)/sinmu;
    double beta_x=one_turn_map[0][1]/sinmu;
    cosmu=0.5*(one_turn_map[2][2]+one_turn_map[3][3]);
    if (fabs(cosmu)>1.) throw std::runtime_error("alpha_y: cosmu larger than zero");
    sinmu=sqrt(1.-cosmu*cosmu);
    double alpha_y=(one_turn_map[2][2]-cosmu)/sinmu;
    double beta_y=one_turn_map[2][3]/sinmu;
    cosmu=0.5*(one_turn_map[4][4]+one_turn_map[5][5]);
    if (fabs(cosmu)>1.) throw std::runtime_error("alpha_z: cosmu larger than zero");
    sinmu=sqrt(1.-cosmu*cosmu);
    double alpha_z=(one_turn_map[4][4]-cosmu)/sinmu;
    double beta_z=one_turn_map[4][5]/sinmu;
   
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
          std::cout<<"beta_x="<<beta_x<<"  alpha_x="<< alpha_x<<std::endl;  
          std::cout<<"beta_y="<<beta_y<<"  alpha_y="<< alpha_y<<std::endl;
          std::cout<<"betaz="<<beta_z<<"  alpha_z="<< alpha_z<<std::endl;
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
        MArray1d clsorbit=stepper_sptr->get_lattice_simulator().get_closed_orbit(0.0);
        
        if (rank==0){
          std::cout<< "chef FracTune x: "<< chef_frac_tunex<< ", EigenTune x: "<< chef_eigen_tunex<<std::endl; 
          std::cout<< "chef FracTune y: "<< chef_frac_tuney<< ", EigenTune y: "<< chef_eigen_tuney<<std::endl;
          
          
          std::cout<< "horizontal chromaticity: "<< horizontal_chromaticity<<std::endl;
          std::cout<< "vertical   chromaticity: "<< vertical_chromaticity<<std::endl;
          std::cout<< "momentum compaction: "<< momentum_compaction<<std::endl;
          std::cout<< "slip factor: "<< slip_factor<<std::endl;
          
          
          std::cout<<"closed orbit: x="<<clsorbit[0]
          <<"  xp="<<clsorbit[1]
          <<"  y="<<clsorbit[2]
          <<"  yp="<<clsorbit[3]
          <<"  ct="<<clsorbit[4]
          <<"  dpop="<<clsorbit[5]
          <<std::endl;
        }
  }// tunes_and_chroms

  stepper_sptr->print_cs_step_betas();
  if (opts.spc_tuneshift) stepper_sptr->cs_step_lattice_functions();

  
  std::vector<double >  actions(6);
  if (!opts.load_bunch) {
        actions= stepper_sptr->get_lattice_simulator().get_stationary_actions(opts.xrms, opts.yrms, opts.zrms/beta);
        if ((actions[0]<=0.) || (actions[1]<=0.) || (actions[2]<=0.)) throw
              std::runtime_error("get_stationary_actions can't satisfy requested moments"); 
        
        if (rank==0) {
          std::cout<<std::endl;
          std::cout<<"actions= ("<< actions[0]<<", "<< actions[1]<<", "<<actions[2]<<")"<<std::endl;
        }
  }

  MArray2d correlation_matrix=get_correlation_matrix(one_turn_map,opts.xrms, opts.yrms, opts.zrms, beta);
  if (rank==0){             
          std::cout<<" correlation matrix="<<multi_array_to_string(correlation_matrix)<<std::endl;
  }
  
  Bunches bunches;
  Commxx_sptr parent_comm_sptr(new Commxx);  
  Commxxs comms(generate_subcomms(parent_comm_sptr, opts.num_bunches));
  for (int i=0;i<opts.num_bunches;++i){
      Commxx_sptr commx=comms[i];
      Bunch_sptr bunch_sptr;
      if (opts.bunch_periodic){
          bunch_sptr=Bunch_sptr(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
          opts.num_macroparticles, opts.num_real_particles, commx,stepper_sptr->get_lattice_simulator().get_bucket_length(),i));   
      }
      else{
        bunch_sptr=Bunch_sptr(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
                            opts.num_macroparticles, opts.num_real_particles, commx));             
      } 
      
      if (commx->has_this_rank()){
          /*if (opts.load_bunch){
            std::stringstream bunch_label;
           // bunch_label<<i;
            bunch_label<<"_m"<<opts.map_order_loaded_bunch<<"_0000";         
            bunch_sptr->read_file("initial_bunch"+bunch_label.str()+".h5");
            if (commx->get_rank()==0) std::cout<<" bunch "<<i<<" loaded from:  "<<"initial_bunch"+bunch_label.str()+".h5"<<std::endl; 
          }
          else{    */         
            
              Random_distribution dist(opts.seed,*commx);
            // populate_6d_stationary_gaussian(dist,  *bunch_sptr, actions, stepper_sptr->get_lattice_simulator());
//               double n_sigma=stepper_sptr->get_lattice_simulator().get_bucket_length()/(6*opts.zrms);
//               if (commx->get_rank()==0) { 
//                   std::cout<<"bunch: "<<i<<": n_sigma for longitudinal gaussian bunch populate: "<< n_sigma<<std::endl;     
//               }      

             
              /*populate_6d_stationary_truncated_longitudinal_gaussian(dist,*bunch_sptr, actions,  n_sigma,
                                                                  stepper_sptr->get_lattice_simulator()); */ 
               MArray1d input_means(boost::extents[6]);
               for(int imean=0;imean<6;++imean){
                 input_means[imean]=0.;
               }
               populate_6d(dist, *bunch_sptr, input_means, correlation_matrix);

              
               if(opts.boxcar) {
                 double boxcar_length=1.*opts.zrms; 
                 populate_longitudinal_boxcar(dist, *bunch_sptr, one_turn_map, boxcar_length);
                 if (commx->get_rank()==0) std::cout<<" longitudinal boxcar distribution"<<std::endl;
               }
               bool kvtransverse=false;
               if( kvtransverse) {
                 double xsize=sqrt(2.)*opts.xrms;
                 double ysize=sqrt(2.)*opts.yrms;
                 double ctrms=opts.zrms/beta;
                 populate_transverseKV_logitudinalGaussian(dist, *bunch_sptr, one_turn_map, xsize, ysize, ctrms);
                 
               }
               
                                                                  
              Lattice_element_slice_sptr first_slice_sptr(*stepper_sptr->get_lattice_simulator().get_slices().begin());         
              if (first_slice_sptr->get_lattice_element().has_string_attribute(
                    "aperture_type")) {
                    std::string   aperture_type = first_slice_sptr->get_lattice_element().get_string_attribute("aperture_type");
                    Aperture_operation_extractor_sptr extractor(
                        stepper_sptr->get_lattice_simulator().get_aperture_operation_extractor_map_sptr()->get_extractor(aperture_type));
                    Aperture_operation_sptr aperture_operation_sptr(extractor->extract(first_slice_sptr));
                    Logger log(0);
                    aperture_operation_sptr->apply(*bunch_sptr,1, log);  
              
              }
//          }   
    }
      
      
      bunches.push_back( bunch_sptr);              
  }   

  Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, stepper_sptr->get_lattice_simulator().get_bucket_length()));
  Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);
 

  for (int i=0;i<bunch_train_sptr->get_size();++i){
    std::stringstream bunch_label;
    bunch_label<<i;
/*
    if (opts.save_bunch){
       bunch_label<<"_m"<<opts.map_order;
       Diagnostics_sptr initial_particles_sptr(new Diagnostics_particles("initial_bunch"+bunch_label.str()+".h5"));
       initial_particles_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
       initial_particles_sptr->update_and_write();     
    }*/
    Diagnostics_sptr diagnostics_full2_sptr(new Diagnostics_full2("step_full2_"+bunch_label.str()+".h5"));
    bunch_train_simulator.add_per_step(i,diagnostics_full2_sptr ); 

//     Diagnostics_sptr diag_operator_sptr(new Diagnostics_full2("operator_diag"+bunch_label.str()+".h5"));
//     bunch_train_simulator.get_diagnostics_actionss().at(i)->add_per_operator(diag_operator_sptr);
  
   //Diagnostics_sptr diagnostics_particles_sptr(new Diagnostics_particles("turn_particles_"+bunch_label.str()+".h5"));
   // bunch_train_simulator.add_per_turn(i,diagnostics_particles_sptr,opts.turn_period);
        
    if(opts.turn_track){
      Diagnostics_sptr diagnostics_bulk_track_sptr(new Diagnostics_bulk_track("bulk_track_"+bunch_label.str()+".h5", 100000));
      bunch_train_simulator.add_per_turn(i,diagnostics_bulk_track_sptr);  

      int grid_z=20;
      int grid_zp=20;
      double z_nsigma=4.0;
      double zp_nsigma=4.0;
      Diagnostics_sptr diagnostics_phase_space_density_sptr(new 
                Diagnostics_phase_space_density("phase_space_density_"+bunch_label.str()+".h5",
                                                grid_z,grid_zp,z_nsigma,zp_nsigma));                                               
      bunch_train_simulator.add_per_turn(i, diagnostics_phase_space_density_sptr);
                                    

    }
    

    if ((bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().has_this_rank())){
//           if (opts.space_charge_be) {
//               if (opts.spc_tuneshift) {
//                 Diagnostics_space_charge_2d_bassetti_erskine_sptr sp_diagnostics_sptr
//                   (new Diagnostics_space_charge_2d_bassetti_erskine("space_charge_be_diagnostics_"+bunch_label.str()+".h5"));
//                 sp_diagnostics_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
//                 spcbe_sptr->add_diagnostics(sp_diagnostics_sptr);
//               }
//           }
          if (opts.space_charge_2dh) { /*space charge 2dh diagnostics to be written*/ }
          if (opts.space_charge_3dh) {
    //         if (opts.spc_tuneshift) {
    //           Diagnostics_space_charge_3d_hockney_sptr sp_diagnostics_sptr
    //                   (new Diagnostics_space_charge_3d_hockney("space_charge_3dh_diagnostics_"+bunch_label.str()+".h5"));
    //               sp_diagnostics_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
    //               spc3dh_sptr->add_diagnostics(sp_diagnostics_sptr);          
    //         }
          }
          if (opts.space_charge_rec) {
                Commxx_sptr comm_spc_sptr=
                    make_optimal_spc_comm(bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm_sptr(),
                            opts.spc_comm_size, opts.equally_spread );
                spcrec_sptr->set_fftw_helper(comm_spc_sptr,opts.equally_spread);
                if (opts.spc_tuneshift) {
                    Diagnostics_space_charge_rectangular_sptr sp_diagnostics_sptr
                        (new Diagnostics_space_charge_rectangular("space_charge_diagnostics_"+bunch_label.str()+".h5"));
                    sp_diagnostics_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
                    spcrec_sptr->add_diagnostics(sp_diagnostics_sptr);
                }
          }
   }//has rank 
   
   
   // adjust means 
    if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().has_this_rank()){
        MArray1d bunch_means=Core_diagnostics::calculate_mean(*bunch_train_simulator.get_bunch_train().get_bunches()[i]); 
        MArray2d_ref particles( bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_local_particles());
        for (int part=0;part<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_local_num();++part){
            particles[part][0] +=-bunch_means[0]+opts.xoffset;
            particles[part][1] +=-bunch_means[1];
            particles[part][2] +=-bunch_means[2]+opts.yoffset;
            particles[part][3] +=-bunch_means[3];
            particles[part][4] +=-bunch_means[4]+beta*opts.zoffset;
            particles[part][5] +=-bunch_means[5];
        }
        bunch_means=Core_diagnostics::calculate_mean(*bunch_train_simulator.get_bunch_train().get_bunches()[i]); 
        MArray1d bunch_stds=Core_diagnostics::calculate_std(*bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means);
        MArray2d bunch_mom2=Core_diagnostics::calculate_mom2(*bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means);
        if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().get_rank()==0) {
            std::cout<<std::endl; 
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
  }//for i
  if ((opts.num_bunches>1) && (rank==0)) std::cout<<"train bunch space="<<bunch_train_simulator.get_bunch_train().get_spacings()[0]<<std::endl;

  
  Propagator propagator(stepper_sptr);
  propagator.set_checkpoint_period(opts.checkpointperiod);
  propagator.set_concurrent_io(opts.concurrentio);
  
  
  
//   for (int i=0;i<bunch_train_sptr->get_size();++i){
//    
//     
//      std::stringstream bunch_label;
//      bunch_label<<i;    
//      Diagnostics_sptr initial_particles_sptr(new Diagnostics_particles("initial_bunch"+bunch_label.str()+".h5"));
//      initial_particles_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
//      initial_particles_sptr->update_and_write();    
//   }
 
  
  

  if ((!opts.longitudinal_modes) && (!opts.transverse_modes)) {
     propagator.propagate(bunch_train_simulator, opts.num_turns, opts.maxturns, opts.verbosity);
  }
  else{
      Excite_mode_actions_sptr exm_actions_sptr;
      if (opts.longitudinal_modes) {
          exm_actions_sptr=Longitudinal_mode_actions_sptr(new 
          Longitudinal_mode_actions(one_turn_map,opts.turn_number_for_action,opts.l_mode_number,0.1));
      }
      else if(opts.transverse_modes){  
        double delta_z=3.;
        exm_actions_sptr=Transverse_mode_actions_sptr(new 
        Transverse_mode_actions(one_turn_map,opts.turn_number_for_action, opts.x_transverse, opts.l_mode_number,
                               // opts.t_mode_number, opts.n_radial, opts.spc_modes,  3., 0.1));
                              opts.t_mode_number, opts.n_radial, opts.spc_modes,opts.mode_from_file, delta_z, opts.delta_t ,opts.mode_file));
                             //opts.t_mode_number, opts.n_radial, opts.spc_modes,opts.mode_from_file,  3., 1 ,opts.mode_file));
      }
      
      propagator.propagate(bunch_train_simulator, *exm_actions_sptr, opts.num_turns, opts.maxturns, opts.verbosity);
    
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
