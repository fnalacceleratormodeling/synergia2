#ifndef MODES_OPTS_H_
#define MODES_OPTS_H_
#include "synergia/utils/logger.h"
#include "synergia/simulation/propagator.h"

class Options
{
  
 
  public:
  

    
    int map_order; 
    int map_order_loaded_bunch;
    bool chef_propagate;
    bool chef_map;
    bool tunes_and_chroms;
    
    bool space_charge_be;
    bool space_charge_2dh;
    bool space_charge_3dh;
    bool space_charge_rec;
    bool if_aperture;
    double aperture;

    double rf_voltage; // "RF cavity voltage in MV"
    int harmon;
  
    double xrms;
    double yrms;
    double zrms;
    
    double xoffset;
    double yoffset;
    double zoffset; 
    
    double k3_e1;
    double k3_e2;
    
    std::vector<int> scgrid; //space_charge grid
       
    int num_steps;
    int steps_per_quad;
    int steps_per_sbend;
    int num_steps_else;
    
    
    int num_bunches;
    bool bunch_periodic;
    bool boxcar;
    int num_macroparticles;
    double num_real_particles;
    long unsigned int seed;
    bool save_bunch;
    bool load_bunch;
    
    bool turn_track;
    int turn_period;
       
    int  spc_comm_size;
    int  equally_spread;
    int  per_host;
    bool spc_tuneshift; // calculates incoherent space charge tune shift
    
    int   num_turns;
    int   checkpointperiod;
    int   maxturns;
    int   concurrentio;
    int   verbosity;
    std::string lattice_file;
    
    std::list<std::string> option_list;
    //std::map<std::string, boost::any> options_map;
   
    int turn_number_for_action;    
    bool longitudinal_modes;
    bool transverse_modes;
    bool x_transverse;
    bool spc_modes;
    bool mode_from_file;
    std::string mode_file;
    int l_mode_number;
    int t_mode_number;
    int n_radial;
    double delta_t;
// default constructors   

    Options(std::string filename="input_options"):
    map_order(1),
    map_order_loaded_bunch(1), 
    space_charge_be(false),
    space_charge_2dh(false),
    space_charge_3dh(false),
    space_charge_rec(false),
    chef_propagate(true),
    chef_map(false),
    tunes_and_chroms(true),
    if_aperture(true),
    aperture(0.06),
    rf_voltage(0.63), 
    harmon(80),  
    xrms(0.008),
    yrms(0.002),
    zrms(1.),  
    xoffset(0.),
    yoffset(0.),
    zoffset(0.),   
    k3_e1(0.),
    k3_e2(0.),
    num_steps(22*100),
    steps_per_quad(1),
    steps_per_sbend(4),
    num_steps_else(22*4),
    num_bunches(2),
    bunch_periodic(1),
    boxcar(false),
    num_macroparticles(1000),
    num_real_particles(5e10),
    seed(13),
    save_bunch(false),
    load_bunch(false), 
    turn_track(0),
    turn_period(1),
    spc_comm_size(32),
    equally_spread(0),
    per_host(0),
    spc_tuneshift(1),
    num_turns(3),
    checkpointperiod(50),
    maxturns(3000),
    concurrentio(8),
    verbosity(1),
    lattice_file("modes_lattice.xml"),
    turn_number_for_action(10),
    longitudinal_modes(0),
    transverse_modes(0),
    x_transverse(1),
    spc_modes(0),    
    mode_from_file(0),
    mode_file("excitation_mode_file.dat"),
    l_mode_number(1),
    t_mode_number(1),
    n_radial(0),
    delta_t(0.3),
    scgrid(3)
    {
      scgrid[0]=32;
      scgrid[1]=32;
      scgrid[2]=64; 
      
    option_list.push_back("map_order");
    option_list.push_back("map_order_loaded_bunch"); 
    option_list.push_back("space_charge_be"); 
    option_list.push_back("space_charge_2dh");
    option_list.push_back("space_charge_3dh");
    option_list.push_back("space_charge_rec");
    option_list.push_back("chef_propagate");
    option_list.push_back("chef_map");
    option_list.push_back("tunes_and_chroms");
    option_list.push_back("if_aperture");
    option_list.push_back("aperture");
    option_list.push_back("rf_voltage"); 
    option_list.push_back("harmon");  
    option_list.push_back("xrms");
    option_list.push_back("yrms");
    option_list.push_back("zrms");  
    option_list.push_back("xoffset");
    option_list.push_back("yoffset");   
    option_list.push_back("zoffset");  
    option_list.push_back("k3_e1");
    option_list.push_back("k3_e2");
    option_list.push_back("num_steps");
    option_list.push_back("steps_per_quad");
    option_list.push_back("steps_per_sbend");
    option_list.push_back("num_steps_else");
    option_list.push_back("num_bunches");
    option_list.push_back("bunch_periodic");
    option_list.push_back("boxcar");
    option_list.push_back("num_macroparticles");
    option_list.push_back("num_real_particles");
    option_list.push_back("seed");
    option_list.push_back("save_bunch");
    option_list.push_back("load_bunch"); 
    option_list.push_back("turn_track");
    option_list.push_back("turn_period");
    option_list.push_back("spc_comm_size");
    option_list.push_back("equally_spread");
    option_list.push_back("per_host");
    option_list.push_back("spc_tuneshift");
    option_list.push_back("num_turns");
    option_list.push_back("checkpointperiod");
    option_list.push_back("maxturns");
    option_list.push_back("concurrentio");
    option_list.push_back("verbosity");
    option_list.push_back("scgrid");
    option_list.push_back("lattice_file");
    option_list.push_back("turn_number_for_action");
    option_list.push_back("longitudinal_modes");
    option_list.push_back("transverse_modes");
    option_list.push_back("x_transverse");
    option_list.push_back("spc_modes");
    option_list.push_back("mode_from_file");
    option_list.push_back("mode_file");
    option_list.push_back("l_mode_number");
    option_list.push_back("t_mode_number");
    option_list.push_back("n_radial");
    option_list.push_back("delta_t");
     
   
    
    
      
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

      std::ifstream rfile(filename.c_str());
      if (rfile.fail()) {
          if (rank==0) std::cout<<"no input options file found, use the default options"<<std::endl;
          return;
      }
      else{
         if (rank==0) std::cout<<"use input options from file: "<<filename<<std::endl;
         std::string line;
          while (!rfile.eof() && rfile.is_open()) {
              getline(rfile,line);
              if ( !line.empty() ){        
                  size_t pos0=line.find_first_not_of(" \t\r\n");
                  if ((pos0 !=std::string::npos) && ((line.at(pos0) != '!') && (line.at(pos0) != '#') )){
                    size_t pos1=line.find_first_of(" =:",pos0);
                    std::string option=line.substr(pos0,pos1-pos0);                  
                    size_t pos2=line.find_first_not_of(" \t\r\n=:",pos1);
                    size_t pos3=line.find_first_of(" ;\t\r\n",pos2);
                    std::string  value= line.substr(pos2,pos3-pos2);  
                    set_option(option,value);
                  }               
                  
              }
          }
      }
     
  }
    





void set_option(std::string option, std::string value){         
    std::list<std::string>::iterator it=find(option_list.begin(), option_list.end(), option);
    if (it==option_list.end()) {
          std::cout<<" option="<<option<<std::endl;
          throw std::runtime_error(" this option is not valid, please modify Options class");
    }
  
     std::stringstream vv;
     if (option=="scgrid")    { 
       size_t pos=value.find_first_of(",");
       std::string v0=value.substr(0,pos);
       size_t  pos1=value.find_first_of(",",pos+1);
       std::string v1=value.substr(pos+1,pos1-pos-1);
       std::string v2(value.begin()+pos1+1,value.end());
       vv.str(v0);
       vv>>scgrid[0];      
       vv.clear();
       vv.str(v1);       
       vv>>scgrid[1];      
       vv.clear();
       vv.str(v2);            
       vv>>scgrid[2]; 
       return;
     }
     
    vv<<value;        
    if (option=="map_order")    vv>>map_order;
    if (option=="map_order_loaded_bunch")    vv>> map_order_loaded_bunch;
    if (option=="space_charge_be")    vv>> space_charge_be;
    if (option=="space_charge_2dh")    vv>> space_charge_2dh;
    if (option=="space_charge_3dh")    vv>> space_charge_3dh;
    if (option=="space_charge_rec")    vv>> space_charge_rec;
    if (option=="chef_propagate")    vv>>chef_propagate;
    if (option=="chef_map")    vv>>chef_map;
    if (option=="tunes_and_chroms")    vv>>tunes_and_chroms;
    if (option=="if_aperture")    vv>>if_aperture;
    if (option=="aperture")    vv>>aperture;
    if (option=="rf_voltage")    vv>> rf_voltage;
    if (option=="harmon")    vv>>  harmon;
    if (option=="xrms")    vv>>xrms;
    if (option=="yrms")    vv>>yrms;
    if (option=="zrms")    vv>> zrms ;
    if (option=="xoffset")    vv>>xoffset;
    if (option=="yoffset")    vv>>yoffset;
    if (option=="zoffset")    vv>> zoffset   ;
    if (option=="k3_e1")    vv>> k3_e1   ;
    if (option=="k3_e2")    vv>> k3_e2   ;
    if (option=="num_steps")    vv>>num_steps;
    if (option=="steps_per_quad")    vv>>steps_per_quad;
    if (option=="steps_per_sbend")    vv>>steps_per_sbend;
    if (option=="num_steps_else")    vv>>num_steps_else;
    if (option=="num_bunches")    vv>>num_bunches;
    if (option=="bunch_periodic")    vv>>bunch_periodic;
    if (option=="boxcar")    vv>>boxcar;
    if (option=="num_macroparticles")    vv>> num_macroparticles;
    if (option=="num_real_particles")    vv>>num_real_particles;
    if (option=="seed")    vv>>seed;
    if (option=="save_bunch")    vv>>save_bunch;
    if (option=="load_bunch")    vv>> load_bunch;
    if (option=="turn_track")    vv>>turn_track;
    if (option=="turn_period")    vv>>turn_period;
    if (option=="spc_comm_size")    vv>>spc_comm_size;
    if (option=="equally_spread")    vv>>equally_spread;
    if (option=="per_host")    vv>>per_host;
    if (option=="spc_tuneshift")    vv>>spc_tuneshift;
    if (option=="num_turns")    vv>>num_turns;
    if (option=="checkpointperiod")    vv>>checkpointperiod;
    if (option=="maxturns")    vv>>maxturns;
    if (option=="concurrentio")    vv>>concurrentio;
    if (option=="verbosity")    vv>>verbosity;  
    if (option=="lattice_file") vv>>lattice_file;
    if (option=="turn_number_for_action") vv>>turn_number_for_action;
    if (option=="longitudinal_modes")vv>>longitudinal_modes;
    if (option=="transverse_modes") vv>>transverse_modes;
    if (option=="x_transverse") vv>>x_transverse;
    if (option=="spc_modes") vv>>spc_modes;
    if (option=="mode_from_file") vv>>mode_from_file;
    if (option=="mode_file") vv>>mode_file;
    if (option=="l_mode_number") vv>>l_mode_number;
    if (option=="t_mode_number") vv>>t_mode_number;
    if (option=="n_radial") vv>>n_radial;
    if (option=="delta_t") vv>>delta_t;
   }
    
   
    void print(Logger  logger=Logger(0,"all_options",false, true) ){
    logger<<"map_order="<<map_order<<std::endl;            
    logger<<"map_order_loaded_bunch="<<map_order_loaded_bunch<<std::endl;    
    logger<<"space_charge_be="<<space_charge_be <<std::endl;
    logger<<"space_charge_2dh="<<space_charge_2dh <<std::endl;
    logger<<"space_charge_3dh="<<space_charge_3dh <<std::endl;
    logger<<"space_charge_rec="<<space_charge_rec <<std::endl;
    logger<<"chef_propagate="<<chef_propagate <<std::endl;
    logger<<"chef_map="<<chef_map<<std::endl;
    logger<<"tunes_and_chroms="<<tunes_and_chroms<<std::endl;
    logger<<"if_aperture="<< if_aperture<<std::endl;
    logger<<"aperture="<< aperture <<std::endl;
    logger<<"rf_voltage="<<rf_voltage <<std::endl; 
    logger<<"harmon="<< harmon<<std::endl;
    logger<<"xrms="<<xrms <<std::endl;
    logger<<"yrms="<< yrms<<std::endl;
    logger<<"zrms="<<zrms <<std::endl;
    logger<<"xoffset="<<xoffset<<std::endl;
    logger<<"yoffset="<<yoffset<<std::endl;
    logger<<"zoffset="<<zoffset <<std::endl;   
    logger<<"k3_e1="<<k3_e1<<std::endl; 
    logger<<"k3_e2="<<k3_e2<<std::endl; 
    logger<<"num_steps="<<num_steps <<std::endl;
    logger<<"steps_per_quad="<<steps_per_quad <<std::endl;
    logger<<"steps_per_sbend="<<steps_per_sbend <<std::endl;
    logger<<"num_steps_else="<<num_steps_else  <<std::endl;
    logger<<"num_bunches="<<num_bunches  <<std::endl;
    logger<<"bunch_periodic="<<bunch_periodic <<std::endl; 
     logger<<"boxcar="<<boxcar<<std::endl; 
    logger<<"num_macroparticles="<<num_macroparticles  <<std::endl;
    logger<<"num_real_particles="<<num_real_particles  <<std::endl;
    logger<<"seed="<<seed  <<std::endl;
    logger<<"save_bunch="<<save_bunch <<std::endl;
    logger<<"load_bunch="<<load_bunch  <<std::endl;
    logger<<"turn_track="<<turn_track <<std::endl;
    logger<<"turn_period="<<turn_period <<std::endl;
    logger<<"spc_comm_size="<<spc_comm_size <<std::endl;
    logger<<"equally_spread="<<equally_spread  <<std::endl;
    logger<<"per_host="<<per_host <<std::endl;
    logger<<"spc_tuneshift="<<spc_tuneshift <<std::endl;
    logger<<"num_turns="<<num_turns  <<std::endl;
    logger<<"checkpointperiod="<<checkpointperiod  <<std::endl;
    logger<<"maxturns="<<maxturns <<std::endl;
    logger<<"concurrentio="<<concurrentio <<std::endl;
    logger<<"verbosity="<<verbosity  <<std::endl;
    logger<<"lattice_file="<<lattice_file<<std::endl;
    logger<<"scgrid="<<scgrid[0]<<","<<scgrid[1]<<","<<scgrid[2]<<std::endl;
    logger<<"turn_number_for_action="<<turn_number_for_action<<std::endl;
    logger<<"longitudinal_modes="<<longitudinal_modes<<std::endl;
    logger<<"transverse_modes="<<transverse_modes<<std::endl;
    logger<<"x_transverse="<<x_transverse<<std::endl;
    logger<<"spc_modes="<<spc_modes<<std::endl;
    logger<<"mode_from_file="<<mode_from_file<<std::endl;
    logger<<"mode_file="<<mode_file<<std::endl;
    logger<<"l_mode_number="<<l_mode_number<<std::endl;
    logger<<"t_mode_number="<<t_mode_number<<std::endl;
    logger<<"n_radial="<<n_radial<<std::endl;
    logger<<"delta_t="<<delta_t<<std::endl;
    }
               
}; 

class Options_resume
{
   public:
     
     bool new_num_turns;
     bool new_maxturns;
     bool new_verbosity;
     
     int num_turns;
     int maxturns;
     int verbosity;
     int checkpointperiod;
     int concurrentio;
     std::string directory;
     
   
     
     std::list<std::string> option_resume_list;
     
     Options_resume(std::string filename="input_resume_options"):  
        new_num_turns(false),
        new_maxturns(false),
        new_verbosity(false),
        num_turns(1),
        maxturns(3000),
        verbosity(1),
        checkpointperiod(50),
        concurrentio(8),
        directory(Propagator::default_checkpoint_dir)
        {
          option_resume_list.push_back("num_turns");
          option_resume_list.push_back("maxturns");
          option_resume_list.push_back("checkpointperiod");
          option_resume_list.push_back("concurrentio");
          option_resume_list.push_back("verbosity");
          option_resume_list.push_back("directory");
          
          
          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
          
          std::ifstream rfile(filename.c_str());
    
          if (rfile.fail()) {
            if (rank==0) std::cout<<"no resume options file found, use the default options"<<std::endl;
            return;
          }
          else{
            if (rank==0) std::cout<<"use resume options from file: "<<filename<<std::endl;
            std::string line;
            while (!rfile.eof() && rfile.is_open()) {
              getline(rfile,line);
              if ( !line.empty() ){        
                size_t pos0=line.find_first_not_of(" \t\r\n");
                if ((pos0 !=std::string::npos) && ((line.at(pos0) != '!') && (line.at(pos0) != '#') )){
                  size_t pos1=line.find_first_of(" =:",pos0);
                  std::string option=line.substr(pos0,pos1-pos0);                  
                  size_t pos2=line.find_first_not_of(" \t\r\n=:",pos1);
                  size_t pos3=line.find_first_of(" ;\t\r\n",pos2);
                  std::string  value= line.substr(pos2,pos3-pos2);  
                  set_option(option,value);
                }               
                
              }
            }
          }
                  
        }
     
     
      void set_option(std::string option, std::string value){         
          std::list<std::string>::iterator it=find(option_resume_list.begin(), option_resume_list.end(), option);
          if (it==option_resume_list.end()) {
                std::cout<<" option="<<option<<std::endl;
                throw std::runtime_error(" this option is not valid, please modify Options_resume class");
          }
        
          std::stringstream vv;          
          vv<<value;        
          if (option=="num_turns"){
                vv>>num_turns;
                new_num_turns=true;
          }
          if (option=="verbosity"){
                vv>>verbosity;
                new_verbosity=true;
          }  
          if (option=="maxturns"){
                vv>>maxturns;
                new_maxturns=true;
          }
          if (option=="checkpointperiod")    vv>>checkpointperiod;
          
          if (option=="concurrentio")    vv>>concurrentio;
          
          if (option=="directory")    vv>>directory;
    
   }
         
     
     void print(Logger  logger=Logger(0,"resume_options",false, true) ){       
        logger<<"num_turns="<<num_turns  <<std::endl;
        logger<<"maxturns="<<maxturns<<std::endl;
        logger<<"checkpointperiod="<<checkpointperiod  <<std::endl;
        logger<<"concurrentio="<<concurrentio <<std::endl;
        logger<<"verbosity="<<verbosity  <<std::endl;
        logger<<"directory="<<directory<<std::endl;
     }
};



#endif /*MODES_OPTS_H_*/

