#include "excite_mode_actions.h"
#include "synergia/utils/multi_array_to_string.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/populate.h"


Excite_mode_actions::Excite_mode_actions()
{
}  

Excite_mode_actions::Excite_mode_actions(int turn_number_for_action):
turn_number_for_action(turn_number_for_action)
{
}  

void
Excite_mode_actions::turn_end_action(Stepper & stepper, Bunch_train & bunch_train, int turn_num)
{
   if (turn_num==this->turn_number_for_action){
      Bunches bunches(bunch_train.get_bunches());
      size_t num_bunches = bunch_train.get_size();
      for (int i = 0; i < num_bunches; ++i)
          if (bunches.at(i)->get_comm().has_this_rank()) {
              if (bunches.at(i)->get_comm().get_rank()==0){
                  std::cout<<" excite_mode applied on bunch "<<i<<"  at end of turn "<<
                  turn_num+1<< " (in fact "<<turn_num<<", counting from 0)"<<std::endl;
              }
                  excite_bunch (*bunches.at(i));                    
          }
     
   }
}

void
Excite_mode_actions::excite_bunch(Bunch &bunch){   
}

template<class Archive>
    void
    Excite_mode_actions::serialize(Archive & ar, const unsigned int version)
    {   
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Propagate_actions);
        ar & BOOST_SERIALIZATION_NVP(turn_number_for_action);      
    }
    

   
template
void
Excite_mode_actions::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Excite_mode_actions::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Excite_mode_actions::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Excite_mode_actions::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
        
        

Excite_mode_actions::~Excite_mode_actions()
{
}  
    
BOOST_CLASS_EXPORT_IMPLEMENT(Excite_mode_actions); 





Longitudinal_mode_actions::Longitudinal_mode_actions()
{
}  

Longitudinal_mode_actions::Longitudinal_mode_actions(MArray2d &one_turn_map, int turn_number_for_action,
                                                   int l_number, double delta_rz ):
Excite_mode_actions(turn_number_for_action),one_turn_map(one_turn_map), l_number(l_number), delta_rz(delta_rz)
{
    double cosmu=0.5*(one_turn_map[4][4]+one_turn_map[5][5]);
    if (fabs(cosmu)>1.) throw std::runtime_error("longitudinal modes: cosmu larger than zero");
    double sinmu=sqrt(1.-cosmu*cosmu);
    alpha_z=(one_turn_map[4][4]-cosmu)/sinmu;
    beta_z=one_turn_map[4][5]/sinmu;
    
}  

void
Longitudinal_mode_actions::excite_bunch(Bunch &bunch){
     if (bunch.get_comm().get_rank()==0){
       std::cout<<"Longitudina one map: alpha="<< alpha_z<<" beta="<<beta_z<<std::endl;
          std::cout<<"Longitudinal_mode_actions:  bunch excited with l mode number "<<l_number
          <<"  and delta_rz="<<delta_rz<<std::endl;
     }
     bunch.convert_to_state(Bunch::fixed_z_lab);
     double z_mean= Core_diagnostics::calculate_z_mean(bunch);
     double zrms=Core_diagnostics::calculate_z_std(bunch,z_mean);
  
     for (int part = 0; part < bunch.get_local_num(); ++part) {
          double cti=bunch.get_local_particles()[part][Bunch::cdt];
          double dppi=bunch.get_local_particles()[part][Bunch::dpop];
          double oci=alpha_z*cti+beta_z*dppi;
          double rri=sqrt(cti*cti+oci*oci); 
          if (rri>1.e-6*zrms){
              double cosmui=cti/rri;
              double sinmui=-oci/rri;
              double rrnew=rri+delta_rz*zrms*cos(l_number*acos(cosmui)); 
              //double rrnew=rri+delta_rz*zrms*sin(l_number*acos(cosmui)); 
              double xnew=rrnew*cosmui;
              double dppnew=(-rrnew*sinmui-alpha_z*xnew)/beta_z;          
              bunch.get_local_particles()[part][Bunch::cdt]=xnew;
              bunch.get_local_particles()[part][Bunch::dpop]=dppnew;
          }    
     }    
}

template<class Archive>
    void
    Longitudinal_mode_actions::serialize(Archive & ar, const unsigned int version)
    {   
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Excite_mode_actions);
        ar & BOOST_SERIALIZATION_NVP(one_turn_map);
        ar & BOOST_SERIALIZATION_NVP(alpha_z); 
        ar & BOOST_SERIALIZATION_NVP(beta_z); 
        ar & BOOST_SERIALIZATION_NVP(delta_rz);
        ar & BOOST_SERIALIZATION_NVP(l_number);      
    }
    

   
template
void
Longitudinal_mode_actions::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Longitudinal_mode_actions::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Longitudinal_mode_actions::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Longitudinal_mode_actions::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
        
        
        



Longitudinal_mode_actions::~Longitudinal_mode_actions()
{
}  
    
BOOST_CLASS_EXPORT_IMPLEMENT(Longitudinal_mode_actions); 

void
Transverse_mode_actions::read_mode()
{
      if (!from_file) return;
      int rank;
      int nlines;
      int nctlines;
      int ndpplines;
      std::vector<double> dpp_coord;
      std::vector<double> ct_coord;
      std::vector<double> dipole_mode;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
      if (rank==0){
        std::cout<<" excitation mode is read from the file: "<<excitation_mode_file_name<<std::endl;
      
        std::ifstream rfile;
        std::string line;
        nlines=0;
        nctlines=0;
        ndpplines=0;
        rfile.open(excitation_mode_file_name.c_str());   
        int nctlines1=0;
        while (!rfile.eof() && rfile.is_open()) {
          getline(rfile,line);
          if ( !line.empty() ){
            size_t pos=line.find_first_not_of(" \t\r\n");
            if (pos !=std::string::npos){
              if (line.at(pos) != '#' ){
                 if (line.at(pos) != '&' ){    
                  nlines++;
                  nctlines1++;
                  std::stringstream ss(line);               
                  double dpp;
                  double ct;
                  double mod;
                  ss>>dpp;
                  ss>>ct;
                  ss>>mod;
                  dpp_coord.push_back(dpp);                                 
                  ct_coord.push_back(ct);
                  dipole_mode.push_back(mod);  
                 }// line not &
                 else{
                   nctlines=nctlines1;
                   nctlines1=0; 
                   ndpplines++;
                 }//line is &
              } // line not #             
            }//pos
          }//line not empty
        }//while
        
        std::cout<<" number of  lines in excitation mode file is="<<nlines<<std::endl;
        std::cout<<" number of ct lines in excitation mode file is="<<nctlines<<std::endl;
        std::cout<<" number of dpp lines in excitation mode file is="<<ndpplines<<std::endl;
        if (nlines != nctlines*ndpplines) 
            throw std::runtime_error(" the mode excitation file is not right");
      }//rank=0
      
       // broadcast dipole_mode, ct_coord, dpp_coord, ndpplines, nctlines, nlines
       
        int error=MPI_Bcast( (void *) &nlines, 1, MPI_INT, 0,  MPI_COMM_WORLD );
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in excite_mode)actions.cc: when broadcasting mode file 1");
        }
        error=MPI_Bcast( (void *) &ndpplines, 1, MPI_INT, 0,  MPI_COMM_WORLD );
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in excite_mode)actions.cc: when broadcasting mode file 2");
        }
        error=MPI_Bcast( (void *) &nctlines, 1, MPI_INT, 0,  MPI_COMM_WORLD );
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in excite_mode)actions.cc: when broadcasting mode file 3");
        }
       
        
        dpp_coord.resize(nlines);
        ct_coord.resize(nlines);
        dipole_mode.resize(nlines);
        error=MPI_Bcast( dpp_coord.data(),  nlines, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in excite_mode)actions.cc: when broadcasting mode file 4");
        }
        error=MPI_Bcast( ct_coord.data(),  nlines, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in excite_mode)actions.cc: when broadcasting mode file 5");
        }
        error=MPI_Bcast(dipole_mode.data(),  nlines, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in excite_mode)actions.cc: when broadcasting mode file 6");
        }
     
      
     
      std::vector<int > grid_shape(3,0);
      grid_shape[0]=ndpplines;
      grid_shape[1]=nctlines;
      grid_shape[2]=1;
      std::vector<double > physical_size(3,0.);
      std::vector<double > physical_offset(3,0.);
      if (ndpplines <2) throw std::runtime_error(" the size of the dpp grid should be larger than 1 ");
      if (nctlines <2) throw std::runtime_error(" the size of the ct grid should be larger than 1 ");
      
      double dpp_cell=(dpp_coord[nlines-1]- dpp_coord[0])/(1.*(ndpplines-1));
      double ct_cell=(ct_coord[nctlines-1] -ct_coord[0])/(1.*(nctlines-1));       
      physical_size[0]= dpp_cell*ndpplines;
      physical_size[1]= ct_cell*nctlines;
      physical_size[2]=0;

      
      rectangular_grid_sptr=Rectangular_grid_sptr(new Rectangular_grid(physical_size,physical_offset,grid_shape,false));
      
      MArray3d_ref mode_values(rectangular_grid_sptr->get_grid_points());
      std::vector<int> shape=rectangular_grid_sptr->get_domain_sptr()->get_grid_shape();
       
      
      for (int i=0;i<shape[0];i++){
        for (int j=0;j<shape[1];j++){
          int index=i* shape[1]+j;
          mode_values[i][j][0]=dipole_mode[index];           
        }
      }         
} 

Transverse_mode_actions::Transverse_mode_actions()
{
}  



Transverse_mode_actions::Transverse_mode_actions(MArray2d &one_turn_map, int turn_number_for_action,
                          bool x_transverse, int l_number, int t_number, int n_radial, bool strong_spc, bool from_file,
                                         double delta_rz, double delta_rt, std::string const & excitation_mode_file_name):
Excite_mode_actions(turn_number_for_action),one_turn_map(one_turn_map), x_transverse(x_transverse),
      l_number(l_number), t_number(t_number), n_radial(n_radial), delta_rz(delta_rz), 
       delta_rt(delta_rt), strong_spc(strong_spc), from_file(from_file), excitation_mode_file_name(excitation_mode_file_name)
{
   
    if (from_file &&   strong_spc) throw std::runtime_error("transverse modes, strong_spc and reading from file cannot be both true");
    double cosmu=0.5*(one_turn_map[4][4]+one_turn_map[5][5]);
    if (fabs(cosmu)>1.) throw std::runtime_error("transverse modes, alpha_z: cosmu larger than zero");
    double sinmu=sqrt(1.-cosmu*cosmu);
    alpha_z=(one_turn_map[4][4]-cosmu)/sinmu;
    beta_z=one_turn_map[4][5]/sinmu;

    if (x_transverse){
      double  cosmu=0.5*(one_turn_map[0][0]+one_turn_map[1][1]); 
      if (fabs(cosmu)>1.) throw std::runtime_error("transverse modes, alpha_x: cosmu larger than zero");
      double sinmu=sqrt(1.-cosmu*cosmu);
      alpha_t=(one_turn_map[0][0]-cosmu)/sinmu;
      beta_t=one_turn_map[0][1]/sinmu;
      
    }
    else{
       double  cosmu=0.5*(one_turn_map[2][2]+one_turn_map[3][3]);        
       if (fabs(cosmu)>1.) throw std::runtime_error("transverse modes, alpha_y: cosmu larger than zero");
       double sinmu=sqrt(1.-cosmu*cosmu);
       alpha_t=(one_turn_map[2][2]-cosmu)/sinmu;
       beta_t=one_turn_map[2][3]/sinmu;
    }
    
    
    if (from_file) read_mode();
      
     
}  

void
Transverse_mode_actions::excite_bunch(Bunch &bunch)
{
    if (bunch.get_comm().get_rank()==0){
          std::cout<<" Transverse one map: alpha_t="<< alpha_t<<" beta_t="<<beta_t
          <<" alpha_z="<< alpha_z<<" beta_z="<<beta_z <<std::endl;      
          if (from_file) {
            std::cout<<" Transverse_mode_actions:  bunch excited with MODE FROM FILE: "<< excitation_mode_file_name
            <<", with delta_rt="<<delta_rt <<",  along x (y) direction "<<x_transverse<<" ("<<!x_transverse<<")"<<std::endl;
          }
          else if(strong_spc){
            std::cout<<" Transverse_mode_actions: STRONG SPACE CHARGE MODES EXCITED with mode="
             <<n_radial <<", with delta_rt="<<delta_rt<<",  along x (y) direction "
             <<x_transverse<<" ("<<!x_transverse<<")"<<std::endl;
          }
          else
          {
            std::cout<<" Transverse_mode_actions:  bunch excited with: l mode number="<<l_number
            <<"  and delta_rz="<<delta_rz<<" n radial number="<<n_radial
            <<", t mode number="<<t_number<<"  and delta_rt="<<delta_rt
            <<",  along x (y) direction "<<x_transverse<<" ("<<!x_transverse<<")"<<std::endl;           
          }
    }
        
     bunch.convert_to_state(Bunch::fixed_z_lab);
     MArray1d means=Core_diagnostics::calculate_mean(bunch); 
     MArray1d stds=Core_diagnostics::calculate_std(bunch,means);
     
     int tt_index, tp_index;
     if (x_transverse){
       tt_index=Bunch::x;
       tp_index=Bunch::xp;      
     }
     else{
       tt_index=Bunch::y;
       tp_index=Bunch::yp;          
     }

     double trms=stds[tt_index];
     double zrms=stds[Bunch::z];
     double mean_z=means[Bunch::z]; 
     
     double dpprms=stds[Bunch::dpop];
     double mean_dpp=means[Bunch::dpop]; 
     
            
     if (from_file){ 
            MArray3d_ref grid_points(rectangular_grid_sptr->get_grid_points());
            std::vector<double > cell_size=rectangular_grid_sptr->get_domain_sptr()->get_cell_size(); 
            std::vector<double > left=rectangular_grid_sptr->get_domain_sptr()->get_left();
            std::vector<int> shape=rectangular_grid_sptr->get_domain_sptr()->get_grid_shape();
            double left_x=left[0];
            double left_y=left[1];
            double cell_x=cell_size[0];
            double cell_y=cell_size[1];
            double shape_x=shape[0];
            double shape_y=shape[1];
       
//            check if the mode is read right from the file           
//           int rank;  
//           MPI_Comm_rank(MPI_COMM_WORLD, &rank);      
//           if (rank==2){               
//               std::cout<<"aa shape dpp="<<shape[0]<<"  shape ct="<<shape[1]<<"  shape 3="<<shape[2] <<std::endl;                   
//               std::cout<<"celle dpp="<<cell_size[0]<<"  cellsize ct="<<cell_size[1]<<"  cellsize 3="<<cell_size[2] <<std::endl;
//               std::vector<double > physical_size=rectangular_grid_sptr->get_domain_sptr()->get_physical_size();
//               std::cout<<"physicsl size[0]="<<physical_size[0]<<"  physical_size ct="<<physical_size[1]<<"  physical_size 3="<<physical_size[2] <<std::endl;
//              
//               std::cout<<" left[0]="<<left[0]<<" left[1]="<<left[1]<<" left[2]="<<left[2]<<std::endl;
//               
//               for (int i=0;i<shape[0];i++){
//                     for (int j=0;j<shape[1];j++){
//                         int index=i* shape[1]+j;
//                         double x,y,z;
//                         rectangular_grid_sptr->get_domain_sptr()->get_cell_coordinates(i,j,0, x,y,z);
//                         std::cout<< x<<"     "<< y <<"     "<< z <<"     "<<grid_points[i][j][0] <<std::endl; 
//                     }
//                     std::cout<<"  &"<<std::endl;
//               }
//           }  
//        
       
       
          for (int part = 0; part < bunch.get_local_num(); ++part) {
                double cti=bunch.get_local_particles()[part][Bunch::cdt];
                double dppi=bunch.get_local_particles()[part][Bunch::dpop];
                double tti=bunch.get_local_particles()[part][tt_index];     
                double x= (dppi-mean_dpp)/dpprms;
                double y=(cti-mean_z)/zrms;         
                double grid_val = interpolate_rectangular_xy(x, y,  left_x, left_y, cell_x, cell_y, shape_x,shape_y, grid_points);
                              
               
                double ttnew=tti-delta_rt*trms*grid_val;
                bunch.get_local_particles()[part][tt_index]=ttnew; 
            
          } 
        
     }
     else{
          for (int part = 0; part < bunch.get_local_num(); ++part) {
                double cti=bunch.get_local_particles()[part][Bunch::cdt];
                double tti=bunch.get_local_particles()[part][tt_index];         
                if (strong_spc)  { 
                      double ttnew=0.;
                      double gaussian_spc_mode=0.;
                      

                      if (n_radial==0) gaussian_spc_mode=1.;// ttnew=tti+delta_rt*trms;

                      if (n_radial==1) {
                        //y = a0*(exp(a1*x)-1)/(exp(a1*x)+1)
                        double  a0=1.61185;
                        double  a1=1.87423;                   
                        gaussian_spc_mode= a0*(exp(a1*(cti-mean_z)/zrms)-1)/(exp(a1*(cti-mean_z)/zrms)+1);                 
                      }
                      if (n_radial==2) {
                        double  a0=2.00444;
                        double  a1=-4.64116;
                        double  a2=-2.38245;
                        double  a3=1.68371;
                        //y=a0+a1/(exp(a2*x-a3)+1)-a1/(exp(a2*x+a3)+1)
                          gaussian_spc_mode=a0+a1/(exp(a2*(cti-mean_z)/zrms-a3)+1.)-a1/(exp(a2*(cti-mean_z)/zrms+a3)+1.);
                          
                                                
                      } 
                        if (n_radial==3) {
                        //y = a0+a1/(exp(a2*x-a3)+1)-a1/(exp(a4*x)+1)+a1/(exp(a2*x+a3)+1)
                            double  a0=2.29763;
                            double a1=-4.59687;
                            double a2=-2.81832;
                            double a3=3.28473;
                            double a4=-4.16831;
                            gaussian_spc_mode=                      
                              a0   +a1/(exp(a2*(cti-mean_z)/zrms-a3)+1.)
                                    -a1/(exp(a4*(cti-mean_z)/zrms)+1.)
                                    +a1/(exp(a2*(cti-mean_z)/zrms+a3)+1.)   ;
                        }
                        if (n_radial==4) {
                        //y= a0+a1/(exp(a2*x-a3)+1)-a1/(exp(a2*x+a3)+1)-a1/(exp(a5*x-a4)+1)+a1/(exp(a5*x+a4)+1)
                          double a0=2.50502;
                          double a1=-5.67112;
                          double a2=2.85366;
                          double a3=3.8121;
                          double a4=1.82341;
                          double a5=4.35194;
                          gaussian_spc_mode=                     
                            a0 + a1/(exp(a2*(cti-mean_z)/zrms-a3)+1)
                                  - a1/(exp(a2*(cti-mean_z)/zrms+a3)+1)
                                  - a1/(exp(a5*(cti-mean_z)/zrms-a4)+1)
                                  + a1/(exp(a5*(cti-mean_z)/zrms+a4)+1)    ;
                        
                        }
                      
                      if (n_radial==-1) {
                        gaussian_spc_mode=-1.55*sin((cti-mean_z)*mconstants::pi/(4.*zrms)); 
                        
                      /*  double dppi=bunch.get_local_particles()[part][Bunch::dpop];
                        //y = a0*(exp(a1*x)-1)/(exp(a1*x)+1)
                        double  a0=1.61185;
                        double  a1=1.87423;                   
                        gaussian_spc_mode= 0.2*a0*(exp(a1*(dppi-mean_dpp)/dpprms)-1)/(exp(a1*(dppi-mean_dpp)/dpprms)+1); */ 
                      }

//                         if (n_radial==-2) {
//                         double dppi=bunch.get_local_particles()[part][Bunch::dpop];
//                         double  a0=2.00444;
//                         double  a1=-4.64116;
//                         double  a2=-2.38245;
//                         double  a3=1.68371;
//                         //y=a0+a1/(exp(a2*x-a3)+1)-a1/(exp(a2*x+a3)+1)
//                         gaussian_spc_mode=a0+a1/(exp(a2*(cti-mean_z)/zrms-a3)+1.)-a1/(exp(a2*(cti-mean_z)/zrms+a3)+1.);
//                         gaussian_spc_mode *=-7.*(dppi-mean_dpp)/dpprms;
//                         }
                      

                      ttnew=tti+delta_rt*trms*gaussian_spc_mode;
                      bunch.get_local_particles()[part][tt_index]=ttnew;               
                }
                else{
                      double dppi=bunch.get_local_particles()[part][Bunch::dpop];
                      double oci=alpha_z*cti+beta_z*dppi;
                      double rri=sqrt(cti*cti+oci*oci); 
                    
              
                      double tpi=bunch.get_local_particles()[part][tp_index];
                      double octi=alpha_t*tti+beta_t*tpi;
                      double qqi=sqrt(tti*tti+octi*octi); 

                      if ((rri>1.e-6*zrms) && (qqi>1.e-6*trms)){
                        double cosmui=cti/rri;
                        double sinmui=-oci/rri;
                        double t_cosmui=tti/qqi;
                        double t_sinmui=-octi/qqi;
      //                   double  qqnew=qqi+delta_rt*trms*cos(t_number*acos(t_cosmui))*delta_rz*zrms
      //                   *cos(n_radial*mconstants::pi*rri/(4.*zrms))*0.5*cos(l_number*acos(cosmui));  
                        double  a0=1.61185;
                        double  a1=1.87423;                   
                        double gaussian_spc_mode= a0*(exp(a1*rri/zrms)-1)/(exp(a1*rri/zrms)+1); 
                        double  qqnew=qqi+delta_rt*trms*cos(l_number*acos(cosmui))*gaussian_spc_mode;
                        
                          double ttnew=qqnew*t_cosmui;
                          double tpnew=(-qqnew*t_sinmui-alpha_t*ttnew)/beta_t;                   
                          bunch.get_local_particles()[part][tt_index]=ttnew;
                          bunch.get_local_particles()[part][tp_index]=tpnew;  
                        
      //                    double ttnew=tti+delta_rt*trms*cos(n_radial*mconstants::pi*rri/(4.*zrms))
      //                    *0.5*(cos(l_number*acos(cosmui))-sin(l_number*acos(cosmui)));
      //                    bunch.get_local_particles()[part][tt_index]=ttnew;
                                        
                    }               
              }    
          }  // for part  
    } //else
//     double jz_clip=5.;
//     clip_gaussian_tails(bunch, one_turn_map, jz_clip);
}


template<class Archive>
void
Transverse_mode_actions::save(Archive & ar, const unsigned int version) const
{           
    ar << BOOST_SERIALIZATION_BASE_OBJECT_NVP(Excite_mode_actions)
        << BOOST_SERIALIZATION_NVP(one_turn_map)
        << BOOST_SERIALIZATION_NVP(alpha_z) 
        << BOOST_SERIALIZATION_NVP(beta_z) 
        << BOOST_SERIALIZATION_NVP(alpha_t) 
        << BOOST_SERIALIZATION_NVP(beta_t)
        << BOOST_SERIALIZATION_NVP(delta_rz)
        << BOOST_SERIALIZATION_NVP(delta_rt)
        << BOOST_SERIALIZATION_NVP(l_number) 
        << BOOST_SERIALIZATION_NVP(n_radial) 
        << BOOST_SERIALIZATION_NVP(t_number)
        << BOOST_SERIALIZATION_NVP(x_transverse)
        << BOOST_SERIALIZATION_NVP(strong_spc)
        << BOOST_SERIALIZATION_NVP(from_file)
        << BOOST_SERIALIZATION_NVP(excitation_mode_file_name);  
       
}


template<class Archive>
void
Transverse_mode_actions::load(Archive & ar, const unsigned int version)
{
            
  ar >> BOOST_SERIALIZATION_BASE_OBJECT_NVP(Excite_mode_actions)
    >> BOOST_SERIALIZATION_NVP(one_turn_map)
    >> BOOST_SERIALIZATION_NVP(alpha_z) 
    >> BOOST_SERIALIZATION_NVP(beta_z) 
    >> BOOST_SERIALIZATION_NVP(alpha_t) 
    >> BOOST_SERIALIZATION_NVP(beta_t)
    >> BOOST_SERIALIZATION_NVP(delta_rz)
    >> BOOST_SERIALIZATION_NVP(delta_rt)
    >> BOOST_SERIALIZATION_NVP(l_number) 
    >> BOOST_SERIALIZATION_NVP(n_radial) 
    >> BOOST_SERIALIZATION_NVP(t_number)
    >> BOOST_SERIALIZATION_NVP(x_transverse)
    >> BOOST_SERIALIZATION_NVP(strong_spc)
    >> BOOST_SERIALIZATION_NVP(from_file)
    >> BOOST_SERIALIZATION_NVP(excitation_mode_file_name);         
    if (from_file) read_mode();        
}


   
template
void
Transverse_mode_actions::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Transverse_mode_actions::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Transverse_mode_actions::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Transverse_mode_actions::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
        
        
        



Transverse_mode_actions::~Transverse_mode_actions()
{
}  
    
BOOST_CLASS_EXPORT_IMPLEMENT(Transverse_mode_actions);





