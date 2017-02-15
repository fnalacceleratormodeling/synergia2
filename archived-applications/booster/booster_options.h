#ifndef BOOSTER_OPTS_H_
#define BOOSTER_OPTS_H_
#include "synergia/utils/logger.h"
#include "synergia/simulation/propagator.h"
#include <boost/any.hpp>



class Options
{
  
 
  public:
    std::string  lattice_file;
    int map_order; 
    int map_order_loaded_bunch;
    bool bpms;
    bool impedance;
    bool space_charge;
    
    bool if_aperture;
    double aperture_f;
    double aperture_d;
    double aperture_l;
    double grid_ref_distance;
    
    double xrms;
    double yrms;
    double zrms;
    
    
    double xoffset;
    double yoffset;
    double zoffset;
  
    std::vector<int> scgrid;    
    std::vector<int> scgrid_l;
    
    int num_steps;
    int steps_per_fmag;
    int steps_per_dmag;
    int num_steps_else;
    
    int num_bunches;
    bool bunch_periodic;
    int num_macroparticles;
    double num_real_particles;
    long unsigned int seed;
    bool load_bunch;
    bool save_bunch;
    bool print_lattice;
    
    
    bool turn_track;
    bool phase_space;
    int turn_period;
    
    
    int  spc_comm_size;
    int  equally_spread;
    int  per_host;
    bool spc_tuneshift;
    
    int   num_turns;
    int   checkpointperiod;
    int   maxturns;
    int   concurrentio;
    int   verbosity;
    
    bool chef_propagate;
    bool chef_map;
    
    
    double rf_voltage;
    int harmon;
    
    bool tunes_and_chroms;
    bool adjust_chromaticity;
    double chrom_h;
    double chrom_v;
    bool adjust_tunes;
    double tune_h;
    double tune_v;
    
    
    std::string wakefile_f;
    std::string wakefile_d;
    std::string  waketype;    
    int  registred_turns;  
    bool full_machine;
    
    std::map<std::string, boost::any> options_map;
    
    
   // std::list<std::string> option_list;
   
    //default constructors   
    Options(std::string filename="input_options"): 
      lattice_file("booster_new.xml"),
      map_order(1),
      map_order_loaded_bunch(1),
      bpms(false),
      impedance(false),
      space_charge(false),
      if_aperture(true),
      aperture_f(0.021),
      aperture_d(0.029),
      aperture_l(0.05715),
      grid_ref_distance(0.01),
      xrms(0.005),
      yrms(0.006),
      zrms(0.4),  
      xoffset(0.),
      yoffset(0.),
      zoffset(0),    
      scgrid(3),
      scgrid_l(3),
      num_steps(24*6),
      steps_per_fmag(1),
      steps_per_dmag(1),
      num_steps_else(24*4),
      num_bunches(1),
      bunch_periodic(1),
      num_macroparticles(1000),
      num_real_particles(5e10),
      seed(13),
      load_bunch(false),
      save_bunch(false), 
      print_lattice(false),
      turn_track(0),      
      phase_space(0),
      turn_period(1),
      spc_comm_size(32),
      equally_spread(0),
      per_host(0),
      spc_tuneshift(1),
      num_turns(1),
      checkpointperiod(50),
      maxturns(3000),
      concurrentio(8),
      verbosity(1),      
      chef_propagate(true),
      chef_map(false),
      rf_voltage(0.6/18.0),
      harmon(84),
      tunes_and_chroms (false),
      adjust_chromaticity(false),
      chrom_h(-15.),
      chrom_v(-10.),
      adjust_tunes (false),
      tune_h(6.755),
      tune_v(6.845),      
      wakefile_f("Fwake.dat"),
      wakefile_d("Dwake.dat"),
      waketype("XLXTYLYTZpp"),      
      registred_turns(15),  
      full_machine(false)
      {
        scgrid[0]=32;     
        scgrid[1]=32;     
        scgrid[2]=64; 
        
        scgrid_l[0]=128;
        scgrid_l[1]=128;
        scgrid_l[2]=64;

          options_map["lattice_file"]= lattice_file;
          options_map["map_order"]=map_order;
          options_map["map_order_loaded_bunch"]=map_order_loaded_bunch;
          options_map["bpms"]=bpms;
          options_map["impedance"]=impedance;
          options_map["space_charge"]=space_charge;
          options_map["if_aperture"]=if_aperture;
          options_map["aperture_f"]=aperture_f;
          options_map["aperture_d"]=aperture_d;
          options_map["aperture_l"]=aperture_l;
          options_map["grid_ref_distance"]=grid_ref_distance;
          options_map["xrms"]=xrms;
          options_map["yrms"]=yrms;
          options_map["zrms"]=zrms;  
          options_map["xoffset"]=xoffset;
          options_map["yoffset"]=yoffset;
          options_map["zoffset"]=zoffset;    
          options_map["scgrid"]=scgrid;
          options_map["scgrid_l"]=scgrid_l;
          options_map["num_steps"]=num_steps;
          options_map["steps_per_fmag"]=steps_per_fmag;
          options_map["steps_per_dmag"]=steps_per_dmag;
          options_map["num_steps_else"]=num_steps_else;
          options_map["num_bunches"]=num_bunches;
          options_map["bunch_periodic"]=bunch_periodic;
          options_map["num_macroparticles"]=num_macroparticles;
          options_map["num_real_particles"]=num_real_particles;
          options_map["seed"]=seed;
          options_map["load_bunch"]=load_bunch;
          options_map["save_bunch"]=save_bunch; 
          options_map["print_lattice"]=print_lattice;
          options_map["turn_track"]=turn_track;
          options_map["phase_space"]= phase_space;
          options_map["turn_period"]=turn_period;
          options_map["spc_comm_size"]=spc_comm_size;
          options_map["equally_spread"]=equally_spread;
          options_map["per_host"]=per_host;
          options_map["spc_tuneshift"]=spc_tuneshift;
          options_map["num_turns"]=num_turns;
          options_map["checkpointperiod"]=checkpointperiod;
          options_map["maxturns"]=maxturns;
          options_map["concurrentio"]=concurrentio;
          options_map["verbosity"]=verbosity;      
          options_map["chef_propagate"]=chef_propagate;
          options_map["chef_map"]=chef_map;
          options_map["rf_voltage"]=rf_voltage;
          options_map["harmon"]=harmon;
          options_map["tunes_and_chroms"]=tunes_and_chroms;
          options_map["adjust_chromaticity"]=adjust_chromaticity;
          options_map["chrom_h"]=chrom_h;
          options_map["chrom_v"]=chrom_v;
          options_map["adjust_tunes"]=adjust_tunes;
          options_map["tune_h"]=tune_h;
          options_map["tune_v"]=tune_v;
          options_map["wakefile_f"]=wakefile_f;
          options_map["wakefile_d"]=wakefile_d;
          options_map["waketype"]=waketype;        
          options_map["registred_turns"]=registred_turns;
          options_map["full_machine"]=full_machine;
              
      
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
        
        

    }// constructor
    

  void set_option(std::string option, std::string value){ 
   if (options_map.find(option) == options_map.end() ){
      std::cout<<" option="<<option<<std::endl;
          throw std::runtime_error(" this option is not valid, please modify Options class");
   }
 

     std::stringstream vv;
     if  (option=="scgrid")    { 
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

    if  (option=="scgrid_l")  { 
       size_t pos=value.find_first_of(",");
       std::string v0=value.substr(0,pos);
       size_t  pos1=value.find_first_of(",",pos+1);
       std::string v1=value.substr(pos+1,pos1-pos-1);
       std::string v2(value.begin()+pos1+1,value.end());
       vv.str(v0);
       vv>>scgrid_l[0];      
       vv.clear();
       vv.str(v1);       
       vv>>scgrid_l[1];      
       vv.clear();
       vv.str(v2);            
       vv>>scgrid_l[2]; 
       return;
     }
     
     
     vv<<value;   
     if (option=="lattice_file") vv>>lattice_file;
     if (option=="map_order")  vv>>map_order;
     if (option=="map_order_loaded_bunch")  vv>>map_order_loaded_bunch;
     if (option=="bpms")  vv>>bpms;
     if (option=="impedance")  vv>>impedance;
     if (option=="space_charge")  vv>>space_charge;
     if (option=="if_aperture")  vv>>if_aperture;
     if (option=="aperture_f")  vv>>aperture_f;
     if (option=="aperture_d")  vv>>aperture_d;
     if (option=="aperture_l")  vv>>aperture_l;
     if (option=="grid_ref_distance")  vv>>grid_ref_distance;
     if (option=="xrms")  vv>>xrms;
     if (option=="yrms")  vv>>yrms;
     if (option=="zrms")  vv>>zrms;  
     if (option=="xoffset")  vv>>xoffset;
     if (option=="yoffset")  vv>>yoffset;
     if (option=="zoffset")  vv>>zoffset;    
     if (option=="num_steps")  vv>>num_steps;
     if (option=="steps_per_fmag")  vv>>steps_per_fmag;
     if (option=="steps_per_dmag")  vv>>steps_per_dmag;
     if (option=="num_steps_else")  vv>>num_steps_else;
     if (option=="num_bunches")  vv>>num_bunches;
     if (option=="bunch_periodic")  vv>>bunch_periodic;
     if (option=="num_macroparticles")  vv>>num_macroparticles;
     if (option=="num_real_particles")  vv>>num_real_particles;
     if (option=="seed")  vv>>seed;
     if (option=="load_bunch")  vv>>load_bunch;
     if (option=="save_bunch")  vv>>save_bunch;
     if (option=="print_lattice")  vv>>print_lattice;
     if (option=="turn_track")  vv>>turn_track;
     if (option=="phase_space")  vv>> phase_space;
     if (option=="turn_period")  vv>>turn_period;
     if (option=="spc_comm_size")  vv>>spc_comm_size;
     if (option=="equally_spread")  vv>>equally_spread;
     if (option=="per_host")  vv>>per_host;
     if (option=="spc_tuneshift")  vv>>spc_tuneshift;
     if (option=="num_turns")  vv>>num_turns;
     if (option=="checkpointperiod")  vv>>checkpointperiod;
     if (option=="maxturns")  vv>>maxturns;
     if (option=="concurrentio")  vv>>concurrentio;
     if (option=="verbosity")  vv>>verbosity;      
     if (option=="chef_propagate")  vv>>chef_propagate;
     if (option=="chef_map")  vv>>chef_map;
     if (option=="rf_voltage")  vv>>rf_voltage;
     if (option=="harmon")  vv>>harmon;
     if (option=="tunes_and_chroms")  vv>>tunes_and_chroms;
     if (option=="adjust_chromaticity") vv>>adjust_chromaticity;
     if (option=="chrom_h") vv>>chrom_h;
     if (option=="chrom_v") vv>>chrom_v;
     if (option=="adjust_tunes") vv>>adjust_tunes;
     if (option=="tune_h") vv>>tune_h;
     if (option=="tune_v") vv>>tune_v;
     if (option=="wakefile_f") vv>>wakefile_f;
     if (option=="wakefile_d") vv>>wakefile_d;
     if (option=="waketype") vv>>waketype;
     if (option=="registred_turns") vv>>registred_turns;
     if (option=="full_machine") vv>> full_machine;
      
   }
    
   
    void print(Logger  logger=Logger(0,"all_options",false, true) ){
     logger<<"lattice_file="<<lattice_file<<std::endl;
     logger<<"map_order="<<map_order<<std::endl;
     logger<<"map_order_loaded_bunch="<<map_order_loaded_bunch<<std::endl;
     logger<<"bpms="<<bpms<<std::endl;
     logger<<"impedance="<<impedance<<std::endl;
     logger<<"space_charge="<<space_charge<<std::endl;
     logger<<"if_aperture="<<if_aperture<<std::endl;
     logger<<"aperture_f="<<aperture_f<<std::endl;
     logger<<"aperture_d="<<aperture_d<<std::endl;
     logger<<"aperture_l="<<aperture_l<<std::endl;
     logger<<"grid_ref_distance="<<grid_ref_distance<<std::endl;
     logger<<"xrms="<<xrms<<std::endl;
     logger<<"yrms="<<yrms<<std::endl;
     logger<<"zrms="<<zrms<<std::endl;  
     logger<<"xoffset="<<xoffset<<std::endl;
     logger<<"yoffset="<<yoffset<<std::endl;
     logger<<"zoffset="<<zoffset<<std::endl;    
     logger<<"num_steps="<<num_steps<<std::endl;
     logger<<"steps_per_fmag="<<steps_per_fmag<<std::endl;
     logger<<"steps_per_dmag="<<steps_per_dmag<<std::endl;
     logger<<"num_steps_else="<<num_steps_else<<std::endl;
     logger<<"num_bunches="<<num_bunches<<std::endl;
     logger<<"bunch_periodic="<<bunch_periodic<<std::endl;
     logger<<"num_macroparticles="<<num_macroparticles<<std::endl;
     logger<<"num_real_particles="<<num_real_particles<<std::endl;
     logger<<"seed="<<seed<<std::endl;
     logger<<"load_bunch="<<load_bunch<<std::endl;
     logger<<"save_bunch="<<save_bunch<<std::endl;   
     logger<<"print_lattice="<<print_lattice<<std::endl;
     logger<<"turn_track="<<turn_track<<std::endl;
     logger<<"phase_space="<<phase_space<<std::endl;
     logger<<"turn_period="<<turn_period<<std::endl;
     logger<<"spc_comm_size="<<spc_comm_size<<std::endl;
     logger<<"equally_spread="<<equally_spread<<std::endl;
     logger<<"per_host="<<per_host<<std::endl;
     logger<<"spc_tuneshift="<<spc_tuneshift<<std::endl;
     logger<<"num_turns="<<num_turns<<std::endl;
     logger<<"checkpointperiod="<<checkpointperiod<<std::endl;
     logger<<"maxturns="<<maxturns<<std::endl;
     logger<<"concurrentio="<<concurrentio<<std::endl;
     logger<<"verbosity="<<verbosity<<std::endl;      
     logger<<"chef_propagate="<<chef_propagate<<std::endl;
     logger<<"chef_map="<<chef_map<<std::endl;
     logger<<"rf_voltage="<<rf_voltage<<std::endl;
     logger<<"harmon="<<harmon<<std::endl;
     logger<<"tunes_and_chroms="<<tunes_and_chroms<<std::endl;
     logger<<"adjust_chromaticity="<< adjust_chromaticity<<std::endl;
     logger<<"chrom_h="<< chrom_h<<std::endl;
     logger<<"chrom_v="<< chrom_v<<std::endl;
     logger<<"adjust_tunes="<< adjust_tunes<<std::endl;
     logger<<"tune_h="<< tune_h<<std::endl;
     logger<<"tune_v="<< tune_v<<std::endl;
     logger<<"wakefile_f="<<wakefile_f<<std::endl;
     logger<<"wakefile_d="<<wakefile_d<<std::endl;
     logger<<"waketype="<<waketype<<std::endl;
     logger<<"registred_turns="<<registred_turns<<std::endl;  
     logger<<"full_machine="<<full_machine<<std::endl; 
     logger<<"scgrid="<<scgrid[0]<<","<<scgrid[1]<<","<<scgrid[2]<<std::endl;
     logger<<"scgrid_l="<<scgrid_l[0]<<","<<scgrid_l[1]<<","<<scgrid_l[2]<<std::endl;      
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

