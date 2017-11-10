#ifndef IOTA_OPTS_H_
#define IOTA_OPTS_H_
#include "synergia/utils/logger.h"
#include "synergia/simulation/propagator.h"
#include <boost/any.hpp>

#endif /*IOTA_OPTS_H_*/

class Options
{
  public:
     std::string  lattice_file;
     int map_order;
     
     bool bpms;
     bool impedance;
     bool space_charge_rec;
     bool space_charge_3dh;
     bool space_charge_2dh;
     
     bool if_aperture;
     bool  coasting_beam;
     double aperture_straight;
     double aperture_bending;
     
     double xrms;
     double yrms;
     double zrms;
     
     
     double xoffset;
     double yoffset;
     double zoffset;
  
     std::vector<int> scgrid_bending;    
     std::vector<int> scgrid_straight; 
     std::vector<int> scgrid_3dh; 
     
     int num_steps;
     int steps_per_sbend;
     int num_steps_straight;
     
     bool turn_track;
     int turn_period;
    
    
     int  spc_comm_size;
     int  equally_spread;
     int  per_host;
     bool spc_tuneshift;
     
     int num_bunches;
     bool bunch_periodic;
     int num_macroparticles;
     double num_real_particles;
     long unsigned int seed;
     bool load_bunch;
     bool save_bunch;
    
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
     
     
     std::string wakefile_straight;
     std::string wakefile_bending;
     std::string  waketype;    
     int  registred_turns;  
     bool full_machine;

     std::map<std::string, boost::any> options_map;



      Options(std::string filename="input_options"): 
      lattice_file("iota_lattice.xml"),
      map_order(1),
      bpms(false),
      impedance(false),
      space_charge_rec(false),
      space_charge_3dh(false),
      space_charge_2dh(false),
      if_aperture(true),
      coasting_beam(false),
      aperture_straight(0.02375),
      aperture_bending(0.025),
      xrms(0.0005),
      yrms(0.0005),
      zrms(0.05), 
      xoffset(0.),
      yoffset(0.),
      zoffset(0),    
      scgrid_bending(3),
      scgrid_straight(3),
      scgrid_3dh(3),
      num_steps(330),
      steps_per_sbend(4),
      num_steps_straight(72*4),
      turn_track(0),
      turn_period(1),
      spc_comm_size(32),
      equally_spread(0),
      per_host(0),
      spc_tuneshift(0),
      num_bunches(1),
      bunch_periodic(1),
      num_macroparticles(1000),
      num_real_particles(1e9),
      seed(13),
      load_bunch(false),
      save_bunch(false),
      num_turns(1),
      checkpointperiod(50),
      maxturns(3000),
      concurrentio(8),
      verbosity(1),
      chef_propagate(true),
      chef_map(false),
      rf_voltage(0),
      harmon(4),
      tunes_and_chroms (false),
      wakefile_straight("IOTA_straight_rw_wake.dat"),
      wakefile_bending("IOTA_bending_rw_wake.dat"),
      waketype("XLYLZ"),      
      registred_turns(15),  
      full_machine(false)
      {
          
          scgrid_bending[0]=64;     
          scgrid_bending[1]=64;     
          scgrid_bending[2]=64; 
        
          scgrid_straight[0]=64;
          scgrid_straight[1]=64;
          scgrid_straight[2]=64; 
          
          scgrid_3dh[0]=64;
          scgrid_3dh[1]=64;
          scgrid_3dh[2]=64; 
          
          options_map["lattice_file"]= lattice_file;
          options_map["map_order"]=map_order;
          options_map["bpms"]=bpms;
          options_map["impedance"]=impedance;
          options_map["space_charge_rec"]=space_charge_rec;
          options_map["space_charge_3dh"]=space_charge_3dh;
          options_map["space_charge_2dh"]=space_charge_2dh;
          options_map["if_aperture"]=if_aperture;
          options_map["coasting_beam"]=coasting_beam;
          options_map["aperture_straight"]=aperture_straight;
          options_map["aperture_bending"]=aperture_bending;
          options_map["xrms"]=xrms;
          options_map["yrms"]=yrms;
          options_map["zrms"]=zrms;
          options_map["xoffset"]=xoffset;
          options_map["yoffset"]=yoffset;
          options_map["zoffset"]=zoffset;    
          options_map["scgrid_bending"]=scgrid_bending;
          options_map["scgrid_straight"]=scgrid_straight;
          options_map["scgrid_3dh"]=scgrid_3dh;
          options_map["num_steps"]=num_steps;
          options_map["steps_per_sbend"]=steps_per_sbend;
          options_map["num_steps_straight"]=num_steps_straight;
          options_map["turn_track"]=turn_track;
          options_map["turn_period"]=turn_period;
          options_map["spc_comm_size"]=spc_comm_size;
          options_map["equally_spread"]=equally_spread;
          options_map["per_host"]=per_host;
          options_map["spc_tuneshift"]=spc_tuneshift;
          options_map["num_bunches"]=num_bunches;
          options_map["bunch_periodic"]=bunch_periodic;
          options_map["num_macroparticles"]=num_macroparticles;
          options_map["num_real_particles"]=num_real_particles;
          options_map["seed"]=seed;
          options_map["load_bunch"]=load_bunch;
          options_map["save_bunch"]=save_bunch; 
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
          options_map["wakefile_straight"]=wakefile_straight;
          options_map["wakefile_bending"]=wakefile_bending;
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
        if  (option=="scgrid_bending")    { 
            size_t pos=value.find_first_of(",");
            std::string v0=value.substr(0,pos);
            size_t  pos1=value.find_first_of(",",pos+1);
            std::string v1=value.substr(pos+1,pos1-pos-1);
            std::string v2(value.begin()+pos1+1,value.end());
            vv.str(v0);
            vv>>scgrid_bending[0];      
            vv.clear();
            vv.str(v1);       
            vv>>scgrid_bending[1];      
            vv.clear();
            vv.str(v2);            
            vv>>scgrid_bending[2]; 
            return;
        }

        if  (option=="scgrid_straight")  { 
            size_t pos=value.find_first_of(",");
            std::string v0=value.substr(0,pos);
            size_t  pos1=value.find_first_of(",",pos+1);
            std::string v1=value.substr(pos+1,pos1-pos-1);
            std::string v2(value.begin()+pos1+1,value.end());
            vv.str(v0);
            vv>>scgrid_straight[0];      
            vv.clear();
            vv.str(v1);       
            vv>>scgrid_straight[1];      
            vv.clear();
            vv.str(v2);            
            vv>>scgrid_straight[2]; 
            return;
        }
        
        
        if  (option=="scgrid_3dh")  { 
            size_t pos=value.find_first_of(",");
            std::string v0=value.substr(0,pos);
            size_t  pos1=value.find_first_of(",",pos+1);
            std::string v1=value.substr(pos+1,pos1-pos-1);
            std::string v2(value.begin()+pos1+1,value.end());
            vv.str(v0);
            vv>>scgrid_3dh[0];      
            vv.clear();
            vv.str(v1);       
            vv>>scgrid_3dh[1];      
            vv.clear();
            vv.str(v2);            
            vv>>scgrid_3dh[2]; 
            return;
        }
        
        
     
        vv<<value;   
        if (option=="lattice_file") vv>>lattice_file;
        if (option=="map_order")  vv>>map_order;
        if (option=="bpms")  vv>>bpms;
        if (option=="impedance")  vv>>impedance;
        if (option=="space_charge_rec")  vv>>space_charge_rec;
        if (option=="space_charge_3dh")  vv>> space_charge_3dh;
        if (option=="space_charge_2dh")  vv>> space_charge_2dh;
        if (option=="if_aperture")  vv>>if_aperture;
        if (option=="coasting_beam")  vv>>  coasting_beam;
        if (option=="aperture_straight")  vv>>aperture_straight;
        if (option=="aperture_bending")  vv>>aperture_bending;
        if (option=="xrms")  vv>>xrms;
        if (option=="yrms")  vv>>yrms;
        if (option=="zrms")  vv>>zrms;  
        if (option=="xoffset")  vv>>xoffset;
        if (option=="yoffset")  vv>>yoffset;
        if (option=="zoffset")  vv>>zoffset;   
        if (option=="num_steps")  vv>>num_steps;
        if (option=="steps_per_sbend")  vv>>steps_per_sbend;
        if (option=="num_steps_straight")  vv>>num_steps_straight;
        if (option=="turn_track")  vv>>turn_track;
        if (option=="turn_period")  vv>>turn_period;
        if (option=="spc_comm_size")  vv>>spc_comm_size;
        if (option=="equally_spread")  vv>>equally_spread;
        if (option=="per_host")  vv>>per_host;
        if (option=="spc_tuneshift")  vv>>spc_tuneshift;
        if (option=="num_bunches")  vv>>num_bunches;
        if (option=="bunch_periodic")  vv>>bunch_periodic;
        if (option=="num_macroparticles")  vv>>num_macroparticles;
        if (option=="num_real_particles")  vv>>num_real_particles;
        if (option=="seed")  vv>>seed;
        if (option=="load_bunch")  vv>>load_bunch;
        if (option=="save_bunch")  vv>>save_bunch;
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
        if (option=="wakefile_straight") vv>>wakefile_straight;
        if (option=="wakefile_bending") vv>>wakefile_bending;
        if (option=="waketype") vv>>waketype;
        if (option=="registred_turns") vv>>registred_turns;
        if (option=="full_machine") vv>> full_machine;
           
     }
     
     
     void print(Logger  logger=Logger(0,"all_options",false, true) ){
         logger<<"lattice_file="<<lattice_file<<std::endl;
         logger<<"map_order="<<map_order<<std::endl;
         logger<<"bpms="<<bpms<<std::endl;
         logger<<"impedance="<<impedance<<std::endl;
         logger<<"space_charge_rec="<<space_charge_rec<<std::endl;
         logger<<"space_charge_3dh="<<space_charge_3dh<<std::endl;
         logger<<"space_charge_2dh="<<space_charge_2dh<<std::endl;
         logger<<"if_aperture="<<if_aperture<<std::endl;
         logger<<"coasting_beam="<<coasting_beam<<std::endl;
         logger<<"aperture_straight="<<aperture_straight<<std::endl;
         logger<<"aperture_bending="<<aperture_bending<<std::endl;
         logger<<"xrms="<<xrms<<std::endl;
         logger<<"yrms="<<yrms<<std::endl;
         logger<<"zrms="<<zrms<<std::endl;
         logger<<"xoffset="<<xoffset<<std::endl;
         logger<<"yoffset="<<yoffset<<std::endl;
         logger<<"zoffset="<<zoffset<<std::endl;
         logger<<"num_steps="<<num_steps<<std::endl;
         logger<<"steps_per_sbend="<<steps_per_sbend<<std::endl;
         logger<<"num_steps_straight="<<num_steps_straight<<std::endl;
         logger<<"turn_track="<<turn_track<<std::endl;
         logger<<"turn_period="<<turn_period<<std::endl;
         logger<<"spc_comm_size="<<spc_comm_size<<std::endl;
         logger<<"equally_spread="<<equally_spread<<std::endl;
         logger<<"per_host="<<per_host<<std::endl;
         logger<<"spc_tuneshift="<<spc_tuneshift<<std::endl;
         logger<<"num_bunches="<<num_bunches<<std::endl;
         logger<<"bunch_periodic="<<bunch_periodic<<std::endl;
         logger<<"num_macroparticles="<<num_macroparticles<<std::endl;
         logger<<"num_real_particles="<<num_real_particles<<std::endl;
         logger<<"seed="<<seed<<std::endl;
         logger<<"load_bunch="<<load_bunch<<std::endl;
         logger<<"save_bunch="<<save_bunch<<std::endl; 
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
         logger<<"wakefile_straight="<<wakefile_straight<<std::endl;
         logger<<"wakefile_bending="<<wakefile_bending<<std::endl;
         logger<<"waketype="<<waketype<<std::endl;
         logger<<"registred_turns="<<registred_turns<<std::endl;  
         logger<<"full_machine="<<full_machine<<std::endl; 
         logger<<"scgrid_bending="<<scgrid_bending[0]<<","<<scgrid_bending[1]<<","<<scgrid_bending[2]<<std::endl;
         logger<<"scgrid_straight="<<scgrid_straight[0]<<","<<scgrid_straight[1]<<","<<scgrid_straight[2]<<std::endl; 
         logger<<"scgrid_3dh="<<scgrid_3dh[0]<<","<<scgrid_3dh[1]<<","<<scgrid_3dh[2]<<std::endl; 
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