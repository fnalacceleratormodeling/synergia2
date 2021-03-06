#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/multi_array_print.h"
#include <boost/archive/archive_exception.hpp>
#include "synergia/utils/multi_array_to_string.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/split_operator_stepper_choice.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/fast_normal_form.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/simulation/populate_stationary.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_phase_space_density.h"
#include "synergia/simulation/diagnostics_normal_form.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/collective/impedance.h"
#include "electronlens_options.h"

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
                std::cerr <<"exception what "<<ex.what()<<std::endl;
                std::cerr << "cxx_electronlens: failed to find "<<opts.lattice_file<<"\n";
                std::cerr << "Run electronlens_xml.py to generate "<<opts.lattice_file<<"\n";
                return;
        }
        if (rank==0) { 
            std::cout<<" lattice file ="<<opts.lattice_file<<std::endl;
        }
        
        
        MPI_Barrier(MPI_COMM_WORLD);   
        double tini = MPI_Wtime();
        
        for (Lattice_elements::const_iterator it =
                lattice_sptr->get_elements().begin();
                it != lattice_sptr->get_elements().end(); ++it) {
            
                std::string element_name=(*it)->get_name();    
                std::string element_type=(*it)->get_type();
       
                (*it)->set_string_attribute("extractor_type", "chef_map");
               
        
                std::string rf_name(element_name.begin(),element_name.begin()+1);    
                if  (rf_name=="r"){
                //std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
                (*it)->set_double_attribute("volt", opts.rf_voltage);
                (*it)->set_double_attribute("lag", 0.);
                (*it)->set_double_attribute("harmon", opts.harmon);
                   //  (*it)->print();
                }
//                 std::string bpm_name(element_name.begin(),element_name.begin()+3);  
//                if (bpm_name=="bpm"){
//                     (*it)->print();
//                }
            
            /*
            if (opts.chef_propagate) {
                (*it)->set_string_attribute("extractor_type", "chef_propagate");
                }
                else if (opts.chef_map) {
                (*it)->set_string_attribute("extractor_type", "chef_map");
                }
                
            if (opts.if_aperture) {
                std::string name1(element_name.begin(),element_name.begin()+1);                      
                if ((element_type=="sbend") && (name1=="m")) {
                    //  (*it)->print(); 
                    //  std::cout<<" elem name1="<<name1<<" type="<<element_type<<std::endl;
                        (*it)->set_string_attribute("aperture_type","rectangular");
                        (*it)->set_double_attribute("rectangular_aperture_width", 2*opts.aperture_bending);
                        (*it)->set_double_attribute("rectangular_aperture_height", 2.*opts.aperture_bending);
                    
                }
                else{
                    (*it)->set_string_attribute("aperture_type","circular");
                    (*it)->set_double_attribute("circular_aperture_radius", opts.aperture_straight);                
                } 
                (*it)->set_string_attribute("aperture_loss","aperture_loss_file");
            }*/

           //   (*it)->print(); 
                
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
        
        
        Space_charge_3d_open_hockney_sptr spc3dh_sptr;    
        
        Dummy_collective_operator_sptr bpm_measure_sptr;
        if (opts.bpms){
            bpm_measure_sptr=Dummy_collective_operator_sptr(new Dummy_collective_operator("bmp_measure"));    
        }
        
        if(opts.space_charge_3dh){
             std::vector<int> grid_shape(opts.scgrid_3dh);
             bool longitudinal_kicks = false;
             bool periodic_z = false;
             double z_period=gamma*lattice_length/opts.harmon;
             bool grid_entire_period = false;
             double nsigma=8.;
             Commxx_divider_sptr commxx_divider_sptr= Commxx_divider_sptr(new Commxx_divider(opts.spc_comm_size, false));
                spc3dh_sptr = Space_charge_3d_open_hockney_sptr( 
                new Space_charge_3d_open_hockney(commxx_divider_sptr, grid_shape,longitudinal_kicks,periodic_z,z_period,grid_entire_period,nsigma));
                spc3dh_sptr->set_green_fn_type(Space_charge_3d_open_hockney::linear);
          
                
             if (rank==0){
                std::cout<<std::endl;
                std::cout<< "3D HOCKNEY SOLVER WITH OPEN BOUNDARIES"<<std::endl;
                std::cout<<"grid for spc 3dh =["<< grid_shape[0]<<", "
                                            << grid_shape[1]<<", "
                                            << grid_shape[2]<<"]"         
                                              <<std::endl;
                std::cout<<"longitudinal kicks ="<<longitudinal_kicks<< std::endl;
                std::cout<<"longitudinal periodicity="<<periodic_z<<std::endl;
                std::cout<<"nsigma="<<nsigma<<std::endl;
                std::cout<<"___________________________________________________________"<<std::endl;                                   
            }
        }
        if ( (opts.space_charge_2dh) || (opts.space_charge_3dh) || (opts.space_charge_rec) || (opts.bpms) ){
                List_choice_map list_choice_map;
        
            Collective_operators operators_quad;
            
            if (opts.space_charge_3dh) operators_quad.push_back(spc3dh_sptr);             
            Kicks  kicks_quad(operators_quad,opts.steps_per_quad);
            
            int steps_per_bpm=1;
            Collective_operators operators_bpm;
            if (opts.bpms) operators_bpm.push_back(bpm_measure_sptr);
            Kicks  kicks_bpm(operators_bpm, steps_per_bpm);                
            
            for (Lattice_elements::const_iterator latt_it =
                lattice_sptr->get_elements().begin();
                latt_it !=lattice_sptr->get_elements().end(); ++latt_it){
                if (opts.space_charge_3dh){
                    if (((*latt_it)->get_name()=="f") || ((*latt_it)->get_name()=="d")) {             
                            // std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
                        list_choice_map[(*latt_it)->get_name()]=kicks_quad;
                    }
                }
                if (opts.bpms) {
                    if((*latt_it)->get_name()=="bpm") {
                        list_choice_map[(*latt_it)->get_name()]=kicks_bpm;
                       // if (rank==0) std::cout<<" BPMS are: "<<(*latt_it)->get_name()<<std::endl;
                    }
                }
                
                                                    
            }
            
            Collective_operators operators_else;               
            if (opts.space_charge_3dh) operators_else.push_back(spc3dh_sptr);
            Kicks  kicks_else(operators_else,opts.num_steps_else);
            list_choice_map["else"]=kicks_else;
          
            stepper_sptr=Stepper_sptr(new Split_operator_stepper_choice(lattice_sptr, opts.map_order, list_choice_map, true)); 
            if (rank==0){
                std::cout<<std::endl; 
                std::cout<<"Split_operator_stepper_choice created"<<std::endl; 
                if (opts.bpms) std::cout<<"measurements at BPMS elements"<<std::endl;
                std::cout<<"steps_per_quad="<<opts.steps_per_quad<<std::endl;
                std::cout<<"num_steps_else="<<opts.num_steps_else<<std::endl;    
                std::cout<<"___________________________________________________________"<<std::endl;
            }                                                   
        }
        else{
            stepper_sptr=Stepper_sptr(new Independent_stepper(lattice_sptr, opts.map_order, opts.num_steps));                 
            if (rank==0) {
                std::cout<<std::endl; 
                std::cout<<"no collective effects, no bpm, independent stepper"<<std::endl;   
                std::cout<<"___________________________________________________________"<<std::endl;
            }
    }
    
      
 
    stepper_sptr->get_lattice_simulator().register_closed_orbit(); 
    stepper_sptr->get_lattice_simulator().set_rf_bucket_length();
   // stepper_sptr->get_lattice_simulator().print_lattice_functions(); 

    if (rank==0) {  
        std::cout<<"stepper:lattice_simulator: map order="<< stepper_sptr->get_lattice_simulator().get_map_order() <<std::endl;
        std::cout<<"stepper: number of steps="<<stepper_sptr->get_steps().size()<<std::endl;
        std::cout<<"stepper rf frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()<<std::endl;
       // stepper_sptr->print();
    }

    stepper_sptr->print_cs_step_betas();
    reference_particle=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_reference_particle();
    lattice_length=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_length();       
    beta = reference_particle.get_beta();
    gamma = reference_particle.get_gamma();
    energy = reference_particle.get_total_energy();
    //  double bunch_sp=stepper_sptr->get_lattice_simulator().get_bucket_length();
    MArray1d clo=stepper_sptr->get_lattice_simulator().get_closed_orbit();   
 
    if (rank==0) { 
        std::cout<<"    ********    stepper  lattice  ************     "<<std::endl;
        std::cout<<std::endl;
        multi_array_print(clo,"closed_orbit");  
//        std::cout<<" bunch spacing="<<bunch_sp<<std::endl;
        std::cout<<" rf bucket length="<<stepper_sptr->get_lattice_simulator().get_bucket_length()<<" ="
         <<stepper_sptr->get_lattice_simulator().get_bucket_length()*1e9/beta/pconstants::c<<" ns"<<std::endl;        
        std::cout<<" rf frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()/1.e6<<" MHz"<<std::endl;     
        std::cout<<" beta="<<beta<<std::endl;
        std::cout<<" gamma="<<gamma<<std::endl;
        std::cout<<" lattice_length="<<lattice_length<<std::endl;
        std::cout<<" closed_orbit_length="<<stepper_sptr->get_lattice_simulator().get_closed_orbit_length()<<std::endl;
        std::cout<<" energy="<<energy<<"GeV"<<std::endl;  
        std::cout<<" reference momentum="<<reference_particle.get_momentum()<<" GeV/c"<<std::endl;        
        std::cout<<std::endl;
        std::cout<<"    ***********************************     "<<std::endl;
        std::cout<<std::endl;
       }

    MArray2d correlation_matrix(boost::extents[6][6]); 
    if (opts.map_order==1){
            MArray2d one_turn_map=stepper_sptr->get_lattice_simulator().get_linear_one_turn_map();
            correlation_matrix=get_correlation_matrix(one_turn_map,opts.xrms, opts.yrms, opts.zrms, beta);
            
            
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
        } 
        double chef_frac_tunex=stepper_sptr->get_lattice_simulator().get_horizontal_tune();
        double chef_frac_tuney=stepper_sptr->get_lattice_simulator().get_vertical_tune();      
        double chef_eigen_tunex=stepper_sptr->get_lattice_simulator().get_horizontal_tune(1);
        double chef_eigen_tuney=stepper_sptr->get_lattice_simulator().get_vertical_tune(1);        
//         double horizontal_chromaticity=stepper_sptr->get_lattice_simulator().get_horizontal_chromaticity();
//         double vertical_chromaticity=stepper_sptr->get_lattice_simulator().get_vertical_chromaticity();                
//         double momentum_compaction=stepper_sptr->get_lattice_simulator().get_momentum_compaction();
//         double slip_factor=stepper_sptr->get_lattice_simulator().get_slip_factor();  
         
        
        if (rank==0){              
            std::cout<< "chef FracTune x: "<< chef_frac_tunex<< ", EigenTune x: "<< chef_eigen_tunex<<std::endl; 
            std::cout<< "chef FracTune y: "<< chef_frac_tuney<< ", EigenTune y: "<< chef_eigen_tuney<<std::endl;                   
           /* std::cout<< "horizontal chromaticity: "<< horizontal_chromaticity<<std::endl;
            std::cout<< "vertical   chromaticity: "<< vertical_chromaticity<<std::endl;
            std::cout<< "momentum compaction: "<< momentum_compaction<<std::endl;
            std::cout<< "slip factor: "<< slip_factor<<std::endl;        */      
            
        } 
        
        Bunches bunches; 
        Commxx_sptr parent_comm_sptr(new Commxx); 
        Commxxs comms(generate_subcomms(parent_comm_sptr, opts.num_bunches));
        for (int i=0;i<opts.num_bunches;++i){
            Commxx_sptr commx=comms[i]; 
            Bunch_sptr bunch_sptr=Bunch_sptr(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
                                                       opts.num_macroparticles, opts.num_real_particles, commx));  
            bunch_sptr->set_bucket_index(i);
            bunch_sptr->set_z_period_length(lattice_sptr->get_length()/opts.harmon) ; 
            if (commx->has_this_rank()){
                if (opts.load_bunch){               
                    bunch_sptr->read_file("el_input_particles.h5");                
                    if (commx->get_rank()==0) std::cout<<" bunch "<<i<<" loaded from:  "<<"el_input_particles.h5"<<std::endl;
                }
                else{
                    if (rank==0) {
                        std::cout<<" input (xrms, yrms, zrms)[m] = ("<< opts.xrms<<", "<< opts.yrms<<", "<< opts.zrms<<")"<<std::endl;                             
                    }
                    Random_distribution dist(opts.seed,*commx);
                    MArray1d input_means(boost::extents[6]);
                    for(int imean=0;imean<6;++imean){
                        input_means[imean]=0.;
                    } 
                    if (opts.map_order==1){
                        populate_6d(dist, *bunch_sptr, input_means, correlation_matrix);
                        
                    }
                    else{ 
                        throw  std::runtime_error("bunch populate for high order maps nor implemented yet");  
                    } 
                } //load bunch            
             
            } //has_this_rank     
            bunches.push_back(bunch_sptr);  
        }// i
        
    
        Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, stepper_sptr->get_lattice_simulator().get_bucket_length()));
        Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

        std::list<int> step_numbers;
        int step_number=0;
        for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it != stepper_sptr->get_steps().end(); ++it) {
            ++step_number;      
            if (opts.bpms){
                for (Operators::const_iterator o_it = (*it)->get_operators().begin(); o_it
                         != (*it)->get_operators().end(); ++o_it) {
                    //  std::cout<<"(*o_it)->get_name()="<< (*o_it)->get_name()<<std::endl; 
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
    
        
        for (int i=0;i<bunch_train_sptr->get_size();++i){
            std::stringstream bunch_label;
            bunch_label<<i;
            

            Diagnostics_sptr diagnostics_full2_sptr(new Diagnostics_full2("step_full2_"+bunch_label.str()+".h5"));
            bunch_train_simulator.add_per_step(i,diagnostics_full2_sptr, step_numbers ); 
            

            Diagnostics_sptr diagnostics_particles_sptr(new Diagnostics_particles("turn_particles_"+bunch_label.str()+".h5"));
            bunch_train_simulator.add_per_turn(i,diagnostics_particles_sptr,opts.turn_period);
            
            
            if(opts.turn_track){
                int nsaved=std::min(10000, opts.num_macroparticles);
                Diagnostics_sptr diagnostics_bulk_track_sptr(new Diagnostics_bulk_track("bulk_track_"+bunch_label.str()+".h5",
                                                                                        nsaved  ));
                bunch_train_simulator.add_per_step(i,diagnostics_bulk_track_sptr, step_numbers );  
        }
            //             if(opts.phase_space){
            //                 int grid_z=20;
            //                 int grid_zp=20;
            //                 double z_nsigma=4.0;
            //                 double zp_nsigma=4.0;
            //                 Diagnostics_sptr diagnostics_phase_space_density_sptr(new 
            //                         Diagnostics_phase_space_density("phase_space_density_"+bunch_label.str()+".h5",
            //                                                         grid_z,grid_zp,z_nsigma,zp_nsigma));                                               
            //                 bunch_train_simulator.add_per_turn(i, diagnostics_phase_space_density_sptr);       
            //             }
            //             
            

            if (opts.space_charge_3dh) {
                if (opts.spc_tuneshift) {
                    Diagnostics_space_charge_3d_hockney_sptr sp_diagnostics_sptr
                        (new Diagnostics_space_charge_3d_hockney("space_charge_3dh_diagnostics_"+bunch_label.str()+".h5"));
                    sp_diagnostics_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);
                    spc3dh_sptr->add_diagnostics(sp_diagnostics_sptr);              
                }           
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
    
  
   
   Propagator propagator(stepper_sptr);
   propagator.set_checkpoint_period(opts.checkpointperiod);
   propagator.set_concurrent_io(opts.concurrentio);
   propagator.propagate(bunch_train_simulator,opts.num_turns, opts.maxturns, opts.verbosity);   
        
}


int
main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   run();
   MPI_Finalize();
    return 0;
}
