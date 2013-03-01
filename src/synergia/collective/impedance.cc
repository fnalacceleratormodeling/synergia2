#include "impedance.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/period.h"
#include "synergia/utils/simple_timer.h"


template<class Archive>
    void
    Bunch_properties::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(x_mean);
        ar & BOOST_SERIALIZATION_NVP(y_mean);
        ar & BOOST_SERIALIZATION_NVP(z_mean);
        ar & BOOST_SERIALIZATION_NVP(realnum);
        ar & BOOST_SERIALIZATION_NVP(bucket_index);
    }

template
void
Bunch_properties::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Bunch_properties::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Bunch_properties::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Bunch_properties::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);










Impedance::Impedance(std::string const & wake_file, std::string const & wake_type, int const  & zgrid,
                    double const & orbit_length, double const & bunchsp, int const nstored_turns,
			                            bool full_machine, std::vector<int > wn):
Collective_operator("impedance"), z_grid(zgrid), nstored_turns(nstored_turns), 
             orbit_length(orbit_length), bunch_spacing(bunchsp), 
              full_machine(full_machine),wn(wn)
{ 
  
  try{
    if (fabs(orbit_length/bunchsp-int(orbit_length/bunchsp))>1e-8)
           throw std::runtime_error("orbit length should divide exacty to bunch_spacing ");
   }
  catch(std::exception const& e){    
        std::cout<<e.what()<<" but the division is "<<orbit_length/bunchsp<< std::endl;   
        MPI_Abort(MPI_COMM_WORLD, 777);
    }  
   
   this->wake_field_sptr=Wake_field_sptr(new Wake_field(wake_file, wake_type)); 
   this->num_buckets=int(orbit_length/bunchsp);   
   construct();        
}  

Impedance::Impedance(std::string const & wake_file, std::string const & wake_type, int const  & zgrid,
                    double const & orbit_length, int const& num_buckets, int const nstored_turns,
			              bool full_machine, std::vector<int > wn):
Collective_operator("impedance"), z_grid(zgrid), nstored_turns(nstored_turns), 
             orbit_length(orbit_length), num_buckets(num_buckets), 
              full_machine(full_machine), wn(wn)
{       
   this->wake_field_sptr=Wake_field_sptr(new Wake_field(wake_file, wake_type)); 
   this->bunch_spacing=orbit_length/num_buckets;
   construct();        
}  

void
Impedance::construct()
{
  try{

    wake_factor=-4.*mconstants::pi*pconstants::rp;   ///N.B. the wakefiled file reads W/(Z_0*L), Z_0=1/(epsilon_0*c)

    if (wn.size() !=3) {
	      wn.resize(3);
	      wn[0]=0;
	      wn[1]=0;
	      wn[2]=0;
       }

    stored_vbunches=std::list< std::vector<Bunch_properties> >();

   // xmom_sptr= boost::shared_ptr<Raw_MArray1d >(new Raw_MArray1d(boost::extents[z_grid]));
    xmom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    ymom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    zdensity_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    xwake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    xwake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    ywake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    ywake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    zwake0_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid])); 
    
    
  }
  catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;   
        MPI_Abort(MPI_COMM_WORLD, 777);
  }  
} 



Impedance::Impedance(Impedance const& impedance)
{

   this->wake_field_sptr=impedance.wake_field_sptr;
   this->z_grid=impedance.z_grid;
   this->nstored_turns=impedance.nstored_turns;
   this->orbit_length=impedance.orbit_length;
   this->num_buckets=impedance.num_buckets;
   this->wake_factor=impedance.wake_factor;
   this->bunch_spacing=impedance.bunch_spacing;
   this->full_machine=impedance.full_machine;
   this->wn=impedance.wn;
   
   // the following data are not copied
   this->stored_vbunches=std::list< std::vector<Bunch_properties> >(); 
   this->xmom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->ymom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->zdensity_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->xwake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->xwake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->ywake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->ywake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   this->zwake0_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
}  


Impedance *
Impedance::clone()
{
    return new Impedance(*this);
}



Impedance::~Impedance(){}

void 
Impedance::set_z_grid(int const  & zgrid)
{
  this->z_grid=zgrid;
}   



Wake_field_sptr 
Impedance::get_wake_field_sptr() const
{
  return wake_field_sptr;
}  


int Impedance::get_z_grid() const { return z_grid;} 
double Impedance::get_orbit_length() const{ return orbit_length;} 
double Impedance::get_wake_factor() const { return wake_factor;} 
double Impedance::get_bunch_spacing() const { return bunch_spacing;}


MArray1d_ref &  Impedance::get_xmom() {return *xmom_sptr;}
MArray1d_ref &  Impedance::get_ymom() {return *ymom_sptr;}
MArray1d_ref &  Impedance::get_zdensity() {return *zdensity_sptr;}
MArray1int_ref & Impedance::get_bin_partition() {return *bin_partition_sptr;}


MArray1d_ref const &  Impedance::get_xmom() const {return *xmom_sptr;}
MArray1d_ref const &  Impedance::get_ymom() const {return *ymom_sptr;}
MArray1d_ref const &  Impedance::get_zdensity() const {return *zdensity_sptr;}
MArray1int_ref const & Impedance::get_bin_partition() const {return *bin_partition_sptr;}

MArray1d_ref &  Impedance::get_xwake_leading() {return *xwake_leading_sptr;}
MArray1d_ref &  Impedance::get_xwake_trailing() {return *xwake_trailing_sptr;}
MArray1d_ref &  Impedance::get_ywake_leading() {return *ywake_leading_sptr;}
MArray1d_ref &  Impedance::get_ywake_trailing() {return *ywake_trailing_sptr;}
MArray1d_ref &  Impedance::get_zwake0() {return *zwake0_sptr;}

MArray1d_ref const &  Impedance::get_xwake_leading() const {return *xwake_leading_sptr;}
MArray1d_ref const &  Impedance::get_xwake_trailing() const {return *xwake_trailing_sptr;}
MArray1d_ref const &  Impedance::get_ywake_leading() const {return *ywake_leading_sptr;}
MArray1d_ref const &  Impedance::get_ywake_trailing() const {return *ywake_trailing_sptr;}
MArray1d_ref const &  Impedance::get_zwake0() const {return *zwake0_sptr;}

bool Impedance::is_full_machine() const { return full_machine;}
int Impedance::get_nstored_turns() const { return nstored_turns;}




void
Impedance::calculate_moments_and_partitions(Bunch & bunch)
{
/// output cell_size_z, xmom, ymom, zdensity
  
    int rank(bunch.get_comm().get_rank());
    
   MArray1d bunchmin(Core_diagnostics::calculate_min(bunch));
   MArray1d bunchmax(Core_diagnostics::calculate_max(bunch));
   // MArray1d bunchmin(Diagnostics::calculate_bunchmin(bunch));
  //  MArray1d bunchmax (Diagnostics::calculate_bunchmax(bunch));
    double z_left= bunchmin[2];
    double z_length=bunchmax[2]-bunchmin[2];       
    cell_size_z= z_length/double(z_grid);
    
    double h = z_length/(z_grid-1.0); 
    if (z_length<= 1.e-14 )   throw
                 std::runtime_error("h and z_length too small ");


    MArray1d_ref xmom(get_xmom());
    MArray1d_ref ymom(get_ymom());
    MArray1d_ref zdensity(get_zdensity());

    int lnum_part=bunch.get_local_num();
    bin_partition_sptr= boost::shared_ptr<MArray1int >(new MArray1int(boost::extents[lnum_part])); 
    MArray1int_ref bin_partition(get_bin_partition());
    
    MArray1d  local_zdensity(boost::extents[z_grid]);
    MArray1d  local_xmom(boost::extents[z_grid]);
    MArray1d  local_ymom(boost::extents[z_grid]);
     
      
    for (int i=0; i<z_grid;  ++i){
        local_zdensity[i]=0.0;
        local_xmom[i]=0.0;
        local_ymom[i]=0.0;
    }



     for (int part = 0;  part < bunch.get_local_num(); ++part) {
         int bin = static_cast<int>((bunch.get_local_particles()[part][4]-z_left)/h);
         if ((bin < z_grid) && (bin >= 0)) {
             local_zdensity[bin] += 1;
             local_xmom[bin] += bunch.get_local_particles()[part][0];
             local_ymom[bin] += bunch.get_local_particles()[part][2];
             bin_partition[part]=bin; //bin_partition(n) is the bin where you find the particle n
         }
         else if ((bin==z_grid) && fabs(bunch.get_local_particles()[part][4]-z_length-z_left)<z_length*1.e-14) { 
             local_zdensity[bin-1] += 1; // put the edge particle in the last bin=z_grid-1
             bin_partition[part]=z_grid-1;
         } 
         else
         {   std::cout << "  z_left  "<<z_left<<"  rank= "<<rank<<std::endl;
          std::cout<<"bunch.get_local_particles()[part][4]="  <<bunch.get_local_particles()[part][4]<<"  rank= "<<rank<<std::endl; 
           std::cout<<"bunch.get_local_particles()[part]0,1,2,3,4,5="  <<bunch.get_local_particles()[part][0]<<
           "  "<<bunch.get_local_particles()[part][1]<<
           "  "<<bunch.get_local_particles()[part][2]<<
           "  "<<bunch.get_local_particles()[part][3]<<
           "  "<<bunch.get_local_particles()[part][4]<<
           "  "<<bunch.get_local_particles()[part][5]<<std::endl;
           
         std::cout<< " particle's id ="<<part<<std::endl; 
         std::cout << "  z_length  "<<z_length<<"  rank= "<<rank<<std::endl;
         std::cout << "(mbs.local_particles(4,n)-z_left)= "<<(bunch.get_local_particles()[part][4]-z_left)<<"  rank= "<<rank<<std::endl;
         std::cout << "bin: " << bin<<"  z_grid="<<z_grid<< "  h=" << h <<"  rank= "<<rank<<std::endl;                
         std::cout << "bunch.get_local_particles()[part][4]-z_length-z_left= "<<fabs(bunch.get_local_particles()[part][4]-z_length-z_left)<<"  rank= "<<rank<<std::endl;
         throw
                 std::runtime_error("particles out of range in calculate_moments_and_partitions ");
         }

    
    }


    int error = MPI_Allreduce(reinterpret_cast<void*>(local_zdensity.origin()), 
			      reinterpret_cast<void*> (zdensity.origin()),
			      z_grid, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance zdensity");
    }
    
    
    MPI_Allreduce(reinterpret_cast<void*>(local_xmom.origin()),
                   reinterpret_cast<void*>(xmom.origin()),
                                           z_grid, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance xmom");
    }
  
    
    MPI_Allreduce(reinterpret_cast<void*>(local_ymom.origin()),
                   reinterpret_cast<void*>(ymom.origin()),
                                           z_grid, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance ymom");
    }    

//     //    dbg is set here 
//    int dbg = 0;    
    for (int k = 0; k < z_grid; ++k) {
       // std::cout<<"zdensity[k]="<<zdensity[k]<<std::endl;
        if (zdensity[k] != 0.0) {
         //   if (dbg) std::cout << "before bin: " << k << " zdensity(k): " << zdensity[k] << " xmom(k): " 
	//	  <<xmom[k] << " ymom(k): " <<  ymom[k] << std::endl;
            xmom[k] /= zdensity[k];
            ymom[k] /= zdensity[k];
	 //    if (dbg) std::cout << "after bin: " << k << " zdensity(k): " << zdensity[k] << " xmom(k): " 
	//	  <<xmom[k] << " ymom(k): " <<  ymom[k] << std::endl;
        } else {
            xmom[k] = 0.0;
            ymom[k] = 0.0;
        }
    }
    
}   









inline int get_zindex_for_wake(double z, double dz, int istart, double zstart)
{ //if  (z< (-istart*istart*dz+zstart)) return -100;
  if (z>=zstart){
    return (static_cast<int>(floor(sqrt((z-zstart)/dz)))) +istart;
  }else
  {
    return (-static_cast<int>(floor(sqrt((zstart-z)/dz)))) +istart;
  } 
}  


void 
Impedance::calculate_kicks()
 {
    
    int zpoints=get_wake_field_sptr()->get_z_coord().size();
    double delta_z=get_wake_field_sptr()->get_delta_z();
    int istart= get_wake_field_sptr()->get_istart();
    double zstart =get_wake_field_sptr()->get_zstart();
    MArray1d_ref z_coord(get_wake_field_sptr()->get_z_coord());
    MArray1d_ref z_wake(get_wake_field_sptr()->get_z_wake());
    MArray1d_ref xw_lead(get_wake_field_sptr()->get_xw_lead());
    MArray1d_ref xw_trail(get_wake_field_sptr()->get_xw_trail());
    MArray1d_ref yw_lead(get_wake_field_sptr()->get_yw_lead());
    MArray1d_ref yw_trail(get_wake_field_sptr()->get_yw_trail());
    

    int registered_turns=stored_vbunches.size();
    int numbunches;
    int num_trains;
    if (registered_turns==0) throw
      std::runtime_error("registered_turns size cannot be zero, probably you propagate a bunch instead of a bunch_train");
    
    numbunches=(*stored_vbunches.begin()).size();
    
    if ((full_machine) && (registered_turns !=0)) {
      num_trains=int(num_buckets/numbunches);
      
      if (fabs(num_buckets/float(numbunches)-num_trains)>1e-8) throw 
	std::runtime_error(
	      "full machine assumes repetitive numer of trains: num_buckets should be divisible to numbunches");
	if (wn[0]<0 || wn[0]>= num_trains ||
	  wn[1]<0 || wn[1]>= num_trains ||
	  wn[2]<0 || wn[2]>= num_trains )
	  throw std::runtime_error(
	      "full machine wave number cannot be smaller than zero or larger than num_trains-1");
    }
    
    
   // std::cout<<" registred turns= "<<registered_turns<<std::endl; 
   // std::cout<<" numbunches= "<<numbunches<<std::endl; 
   // std::cout<<" num_buckets= "<<num_buckets<<std::endl;
   //  std::cout<<" num_trains= "<<num_trains<<std::endl;
   //  std::cout<<"wn="<<wn[0]<<","<<wn[1]<<","<<wn[2]<<std::endl;
   
   MArray1d_ref const xmom(get_xmom());
   MArray1d_ref const ymom(get_ymom());
   MArray1d_ref const zdensity(get_zdensity());
   MArray1d_ref xwake_leading(get_xwake_leading());
   MArray1d_ref xwake_trailing(get_xwake_trailing());
   MArray1d_ref ywake_leading(get_ywake_leading());
   MArray1d_ref ywake_trailing(get_ywake_trailing()); 
   MArray1d_ref zwake0(get_zwake0());
   
   
    for (int i = 0; i < z_grid; ++i){
        xwake_leading[i]=0.; 
	xwake_trailing[i]=0. ;
        ywake_leading[i] =0.; 
	ywake_trailing[i]=0. ;
	zwake0[i] =0.;	
      // in-bunch impedance 
        for (int j = 0; j < z_grid; ++j){
       // for (int j = i+1; j < z_grid; ++j){
            double zji=(j-i)*cell_size_z;
            if (zji>=z_coord[0]) {// at small distance no impedance is considered
                // below it is assumed the wake function is stored using a quadratic grid           
	        int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                double xwl(0.), xwt(0.), ywl(0.), ywt(0.), zw(0.);
                    if (iz+1 < zpoints) {
                    xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
		    xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
		    ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
		    ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                    zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]); 
                    }
                xwake_leading[i]  +=zdensity[j]*N_factor*xmom[j]*xwl; 
		xwake_trailing[i]  += zdensity[j]*N_factor*xwt;
                ywake_leading[i]  += zdensity[j]*N_factor*ymom[j]*ywl; 
		ywake_trailing[i]  += zdensity[j]*N_factor*ywt;
	        zwake0[i] += zdensity[j]*N_factor*zw;				
            }          
         }
        
        
        std::list< std::vector<Bunch_properties> >::const_iterator it=stored_vbunches.begin(); // stored_vbunches.begin() stores the bunches info at 
	                                                                                  // at the moment
        /// bucket 0 is in front of bucket 1, which is in front of bucket 2, etc...
        double z_to_edge=(z_grid-i-1)*cell_size_z;
        for (int ibunch= 0; ibunch<numbunches; ++ibunch){            
	  //  double xwl(0.), xwt(0.), ywl(0.), ywt(0.), zw(0.);
            int ibucket=(*it)[ibunch].bucket_index;
            if(ibucket<bunch_bucket) {///  same turn, the leading buckets effect    
            double  zji=z_to_edge+bunch_spacing*(bunch_bucket-ibucket) +((*it)[ibunch].z_mean-bunch_z_mean);
          //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));  
	    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
            if ((iz+1 < zpoints) && (iz>0)) {
		    double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
		    double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
		    double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
		    double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                    double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]); 
		    
		    xwake_leading[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
		    xwake_trailing[i]  += (*it)[ibunch].realnum*xwt;
		    ywake_leading[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
		    ywake_trailing[i]  +=(*it)[ibunch].realnum*ywt;
		    zwake0[i] += (*it)[ibunch].realnum*zw;
                }
            }
        } // ibunch loop

         if (full_machine) { // it assumes that all the other trains are in front 
 	   for (int itrain= 1; itrain<num_trains; ++itrain){
 	     double wnx=cos(2.*mconstants::pi*wn[0]*itrain/double(num_trains));
 	     double wny=cos(2.*mconstants::pi*wn[1]*itrain/double(num_trains));
 	     double wnz=cos(2.*mconstants::pi*wn[2]*itrain/double(num_trains));
 	     for (int ibunch= 0; ibunch<numbunches; ++ibunch){		  
 		   double  zji=z_to_edge+bunch_spacing*numbunches*itrain+bunch_spacing*(bunch_bucket-ibunch)+
 		    ((*it)[ibunch].z_mean*wnz-bunch_z_mean);
		    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0])))); 
 		   int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
 		   if ((iz+1 < zpoints) && (iz>0)) {
			  double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
			  double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);    
			  
			  xwake_leading[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
			  xwake_trailing[i]  += (*it)[ibunch].realnum*xwt;
			  ywake_leading[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
			  ywake_trailing[i]  +=(*it)[ibunch].realnum*ywt;
			  zwake0[i] += (*it)[ibunch].realnum*zw;
                   }
 	      }
 	   }
	}  // full_machine
        
         if (registered_turns>1) {
            ++it; ///previous turn, following buckets effect
            for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                 int ibucket=(*it)[ibunch].bucket_index;
                 if(ibucket>=bunch_bucket) {///  following buckets effect    
                    double  zji=z_to_edge+bunch_spacing*(bunch_bucket-ibucket)+orbit_length +((*it)[ibunch].z_mean-bunch_z_mean);
                  //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0])))); 
		    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
			  double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
			  double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);    
			  
			  xwake_leading[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
			  xwake_trailing[i]  += (*it)[ibunch].realnum*xwt;
			  ywake_leading[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
			  ywake_trailing[i]  +=(*it)[ibunch].realnum*ywt;
			  zwake0[i] += (*it)[ibunch].realnum*zw;
			  
                     }
                 }
            } // ibunch loop	 
	}//registered_turns>1
         
    } // i loop
    /// it is not necessary to have a loop over i at larger distances, since the effect is negligible
    if (registered_turns>1) {
        double xwake_l=0.;
        double xwake_t=0.;
	double ywake_l=0.;
        double ywake_t=0.;
        double zwake_0=0.; 
              
        
        std::list< std::vector<Bunch_properties> >::const_iterator it;
	std::list< std::vector<Bunch_properties> >::const_iterator jt=stored_vbunches.begin();
	++jt;       
	int iturn;
        for (it=jt,  iturn=1; it !=stored_vbunches.end(); ++it, ++iturn){
            for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                int ibucket=(*it)[ibunch].bucket_index;
               // if(((ibucket<bunch_bucket) && (it==jt))///  finishing the previous turn, for the buckets ahead
               //                 ||  (it!=jt))  /// previous turns effects
                if(((ibucket<bunch_bucket) && (iturn==1))///  finishing the previous turn, for the buckets ahead
                                ||  (iturn>1))  /// previous turns effects                
                                
                {
                    double  zji=bunch_spacing*(bunch_bucket-ibucket)+orbit_length*iturn +((*it)[ibunch].z_mean-bunch_z_mean);
                  //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));  
		    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
			  double  xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double  xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double  ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double  ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
			  double  zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);
			  
		          xwake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl;
			  xwake_t += (*it)[ibunch].realnum*xwt;
			  ywake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
			  ywake_t += (*it)[ibunch].realnum*ywt; 
			  zwake_0 += (*it)[ibunch].realnum*zw;
                    }
                } 
            }
            
            if (full_machine) {
	      for (int itrain= 1; itrain<num_trains; ++itrain){
		double wnx=cos(2.*mconstants::pi*wn[0]*itrain/double(num_trains));
		double wny=cos(2.*mconstants::pi*wn[1]*itrain/double(num_trains));
		double wnz=cos(2.*mconstants::pi*wn[2]*itrain/double(num_trains));
		for (int ibunch= 0; ibunch<numbunches; ++ibunch){
		     double  zji=orbit_length*iturn+bunch_spacing*numbunches*itrain+bunch_spacing*(bunch_bucket-ibunch)+
 		     ((*it)[ibunch].z_mean*wnz-bunch_z_mean);
		  //    int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));
		      int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
		      if ((iz+1 < zpoints) && (iz>0)) {
			  double  xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double  xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double  ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
			  double  ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
			  double  zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);
			  
		          xwake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl;
			  xwake_t += (*it)[ibunch].realnum*xwt;
			  ywake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
			  ywake_t += (*it)[ibunch].realnum*ywt; 
			  zwake_0 += (*it)[ibunch].realnum*zw;
                    }		      
		}//itrain
	      } //ibunch
	    }  // full_machine
      } //iturn
      for (int i = 0; i < z_grid; ++i){
	xwake_leading[i] +=  xwake_l;
	xwake_trailing[i] += xwake_t;
	ywake_leading[i] +=  ywake_l; 
	ywake_trailing[i] += ywake_t ;
	zwake0[i] += zwake_0;    
      }         
   }//registred_turns>1
}


void 
Impedance::apply_impedance_kick(Bunch & bunch, double wake_factor)
{
  
 MArray1int_ref const bin_partition(get_bin_partition());
 MArray1d_ref const  xwake_leading(get_xwake_leading());
 MArray1d_ref const xwake_trailing(get_xwake_trailing());
 MArray1d_ref const  ywake_leading(get_ywake_leading());
 MArray1d_ref const  ywake_trailing(get_ywake_trailing()); 
 MArray1d_ref const zwake0(get_zwake0());
 
 
 for (int part = 0; part < bunch.get_local_num(); ++part) {
        double xkick=0., ykick=0., zkick=0.;
        int bin=bin_partition[part];  // bin_partition(n) is the bin where you find the particle n 
//            if ((bin>=z_grid) || (bin<0))  { std::cout<<"bin="<<bin<<"z_grid="<<z_grid<<std::endl;
//         throw
//         std::runtime_error("something wrong with bining");} 

        xkick=xwake_leading[bin]+xwake_trailing[bin]*bunch.get_local_particles()[part][0];	
        ykick=ywake_leading[bin]+ywake_trailing[bin]*bunch.get_local_particles()[part][2];
	zkick = zwake0[bin];
    
       
        bunch.get_local_particles()[part][1] += wake_factor*xkick;   
        bunch.get_local_particles()[part][3]  += wake_factor*ykick;
        bunch.get_local_particles()[part][5]  += wake_factor*zkick;
      
      
    }
}








void
Impedance::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{
   double t;
   t = simple_timer_current();
   bunch.convert_to_state(Bunch::fixed_t_lab);  
   calculate_moments_and_partitions(bunch);      
      
//       std::ofstream file;
//       file.open("zdensity.dat");
// 	for (int i = 0; i < z_grid; ++i){
// 	file<<i<<"   "<<(*zdensity_sptr)[i]<<"   "<<(*xmom_sptr)[i]<<"   "<<(*ymom_sptr)[i]<<std::endl;
// 	}  
//       file.close();
//       abort(); 
   N_factor=bunch.get_real_num()/bunch.get_total_num();
   bunch_z_mean=Core_diagnostics::calculate_z_mean(bunch);
   bunch_bucket=bunch.get_bucket_index(); 
   calculate_kicks();

   double gamma = bunch.get_reference_particle().get_gamma();
   double beta= bunch.get_reference_particle().get_beta();
   double w_f=get_wake_factor()*time_step/(gamma*beta);  
      
   apply_impedance_kick(bunch,  w_f);
   t = simple_timer_show(t, "impedance apply");
   
}



void
Impedance::store_bunches_data(Bunch_train & bunch_train)
{       
  
    Bunches bunches(bunch_train.get_bunches());
    size_t num_bunches = bunch_train.get_size();   
    Bunch_properties bi;
    std::vector<Bunch_properties> vbi_local(0);
    std::vector<Bunch_properties> vbi(num_bunches);    
    for (int i = 0; i < num_bunches; ++i)
        if (bunches.at(i)->get_comm().has_this_rank()) {
	    Bunch_sptr bunch_sptr=bunches.at(i);
            bunch_sptr->convert_to_state(Bunch::fixed_t_lab);
	    MArray1d bunch_means=Core_diagnostics::calculate_mean(*bunch_sptr);
	    if (full_machine)  
	      if  (bunch_sptr->get_bucket_index() != i) 				
	                   throw std::runtime_error("for full_machine the buckets have to be occupied in order");		     		    
            bi.x_mean=bunch_means[0];
	    bi.y_mean=bunch_means[2];
	    bi.z_mean=bunch_means[4];
	    bi.realnum=bunch_sptr->get_real_num();
	    bi.bucket_index=bunch_sptr->get_bucket_index();  	
	    if  (bunch_sptr->get_comm().get_rank()==0)   vbi_local.push_back(bi);
                     ///only the rank 0 of every communicator sends the bi to all others
        }  
  
	    
	     
           
	MPI_Datatype Bunch_properties_type;
	MPI_Aint lb, extent;
	MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent); 
	MPI_Datatype type[2] = {MPI_DOUBLE, MPI_INT};
	int blocklen[2] = {4,1};
	MPI_Aint disp[2];
	disp[0]=0;
	disp[1]=4*extent;
	MPI_Type_create_struct(2,blocklen, disp, type, &Bunch_properties_type);
	MPI_Type_commit(&Bunch_properties_type); 
                                     
        int size_parent_comm=bunch_train.get_parent_comm_sptr()->get_size();	
	std::vector<int > counts(bunch_train.get_proc_counts_for_impedance());
	std::vector<int > offsets(bunch_train.get_proc_offsets_for_impedance());
	
	
	int error = MPI_Allgatherv(reinterpret_cast<void*>(&vbi_local[0]), vbi_local.size(), Bunch_properties_type,  
					reinterpret_cast<void*>(&vbi[0]), &counts[0], &offsets[0], 
					Bunch_properties_type, MPI_COMM_WORLD);
	if (error != MPI_SUCCESS) {
	  throw std::runtime_error("Impedance::store_bunches_data: MPI error in MPI_Allgatherv");
	} 
		    
	MPI_Type_free(&Bunch_properties_type);
	stored_vbunches.push_front(vbi);          
	if (stored_vbunches.size()>nstored_turns) stored_vbunches.pop_back();
         

}


void
Impedance::apply(Bunch_train & bunch_train, double time_step, Step & step,
        int verbosity, Train_diagnosticss const& per_operation_diagnosticss,
        Logger & logger)
{ 
    store_bunches_data(bunch_train);
    Bunches bunches(bunch_train.get_bunches());
    size_t num_bunches = bunch_train.get_size();
    for (int i = 0; i < num_bunches; ++i)
        if (bunches.at(i)->get_comm().has_this_rank()) {
            apply(*bunches.at(i), time_step, step, verbosity,logger);
        }  
} 
   
    
