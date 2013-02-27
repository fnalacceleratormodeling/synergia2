#include "wake_field.h"
#include <fstream>
#include <stdexcept>
#include <sstream>

Wake_field::Wake_field(std::string const & wake_file, std::string const & wake_type):
wake_file(wake_file), wake_type(wake_type)
{
  
  try{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    int size_wake;
    if (rank==0){
      // read the wakes from the file wake_file     
      // for parallel plates geometry wake file should be written as a four column file such, containing wakes functions such:
      //  z[m]        Wz_trx/Z_0/L[1/(m^2*s]       Wz_try/Z_0/L[1/(m^2*s]        Wz_l/Z_0/L[1/(ms)]   
      // the lines starting with "#" in the file are skipped and can be  used for comments 
	    std::ifstream rfile;
	    std::string line;
	    rfile.open(wake_file.c_str());
	    int num_columns_prev(-1);
	    while (!rfile.eof() && rfile.is_open()) {
	      std::vector<double> temp_wake;
	      getline(rfile,line);
	      
	      if ( !line.empty() ){  
		  int pos=line.find_first_not_of(" \t\r\n");
		  if (pos !=std::string::npos){
		      if (line.at(pos) != '#' ){
			    std::stringstream ss(line);
			    double column;
			    int num_columns(0);
			    while (ss>>column){
			      temp_wake.push_back(column);
			      num_columns++;			  
			    }
			    
			    if (!((num_columns_prev ==-1) || (num_columns_prev == num_columns)))  throw
				std::runtime_error(" the number of columns in the wake file is not the same on all lines");
			    num_columns_prev=num_columns;
			    if (num_columns==1) throw
				  std::runtime_error(" the wake file has only one column");
			    
			    if (num_columns==2) {
				  if (temp_wake.size()!=2) throw std::runtime_error(" temp_wake size should be 2 in this case");
				  z_coord.push_back(temp_wake[0]);				 				 
				  if (get_wake_type()=="Z") {
				      z_wake.push_back(temp_wake[1]);
				  }  				 
				  /* else if(get_wake_type()=="XL") {
				      xw_lead.push_back(temp_wake[1]);
				  }
				  else if(get_wake_type()=="XT") {
				      xw_trail.push_back(temp_wake[1]);
				  }
				    else if(get_wake_type()=="YL") {
				      yw_lead.push_back(temp_wake[1]);
				  }
				  else if(get_wake_type()=="YT") {
				      yw_trail.push_back(temp_wake[1]);
				  } */
				  else if (get_wake_type()=="XLYL") {				   
				      xw_lead.push_back(temp_wake[1]);
				      yw_lead.push_back(temp_wake[1]);
				  } 
				  else{
				    throw
				      std::runtime_error("invalid specification of the wake type for 2 columns wake file");
				  } 
				  continue;
			    }
			    
			    if (num_columns==3) {
				  if (temp_wake.size()!=3) throw std::runtime_error(" temp_wake size should be 3 in this case");
				  z_coord.push_back(temp_wake[0]);				 				 
				  /*if (get_wake_type()=="XLZ") {
				      xw_lead.push_back(temp_wake[1]);
				      z_wake.push_back(temp_wake[2]);				    
				  }  
				  if (get_wake_type()=="XTZ") {
				      xw_trail.push_back(temp_wake[1]);
				      z_wake.push_back(temp_wake[2]);				    
				  }  
				  if (get_wake_type()=="YLZ") {
				      yw_lead.push_back(temp_wake[1]);
				      z_wake.push_back(temp_wake[2]);				    
				  } 
				  if (get_wake_type()=="YTZ") {
				      yw_trail.push_back(temp_wake[1]);
				      z_wake.push_back(temp_wake[2]);				    
				  } 
				  if (get_wake_type()=="XLYL") {
				      xw_lead.push_back(temp_wake[1]);
				      yw_lead.push_back(temp_wake[2]);					    
				  } 
				  if (get_wake_type()=="XTYT") {
				      xw_trail.push_back(temp_wake[1]);
				      yw_trail.push_back(temp_wake[2]);					    
				  } 
				  if (get_wake_type()=="XLXT") {
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(temp_wake[2]);					    
				  } 
				  if (get_wake_type()=="YLYT") {
				      yw_lead.push_back(temp_wake[1]);
				      yw_trail.push_back(temp_wake[2]);					    
				  } 			*/	 
				  if (get_wake_type()=="XLYLZ") {
				      xw_lead.push_back(temp_wake[1]);
				      yw_lead.push_back(temp_wake[1]);
				      z_wake.push_back(temp_wake[2]);
				  } 
				  else if (get_wake_type()=="XLXTYLYT") {
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(temp_wake[2]);
				      yw_lead.push_back(temp_wake[1]);
				      yw_trail.push_back(temp_wake[2]);				   
				  } 	
				  else{
				    throw
				      std::runtime_error("invalid specification of the wake type for 3 columns wake file");
				  } 
				  continue;
			    }
			    
			    
			    if (num_columns==4) {
				  if (temp_wake.size()!=4) throw std::runtime_error(" temp_wake size should be 4 in this case");
				  z_coord.push_back(temp_wake[0]);				 				 
				  
				  if (get_wake_type()=="XLXTYLYTZ") {
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(temp_wake[2]);
				      yw_lead.push_back(temp_wake[1]);
				      yw_trail.push_back(temp_wake[2]);
				      z_wake.push_back(temp_wake[3]);
				  } 
				  else if (get_wake_type()=="XLYLZ") {
				    xw_lead.push_back(temp_wake[1]);
				    yw_lead.push_back(temp_wake[2]);
				    z_wake.push_back(temp_wake[3]);
				  } 
				  else  if (get_wake_type()=="XLXTYLYTZpp") { // this is for old wake files for parallel planes geometry
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(-temp_wake[1]);
				      yw_lead.push_back(temp_wake[2]);
				      yw_trail.push_back(temp_wake[1]);
				      z_wake.push_back(temp_wake[3]);
				  } 
				  else{
				    throw
				      std::runtime_error("invalid specification of the wake type for 4 columns wake file");
				  } 
				  continue;
			    }
			    
			    if (num_columns==5) {
				  if (temp_wake.size()!=5) throw std::runtime_error(" temp_wake size should be 5 in this case");
				  z_coord.push_back(temp_wake[0]);				 				 
				  
				  if (get_wake_type()=="XLXTYLYT") {
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(temp_wake[2]);
				      yw_lead.push_back(temp_wake[3]);
				      yw_trail.push_back(temp_wake[4]);				    
				  } 
				    else  if (get_wake_type()=="XLXTYLYTZpp"){
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(-temp_wake[1]);  
				      yw_lead.push_back(temp_wake[2]);
				      yw_trail.push_back(temp_wake[3]);
				      z_wake.push_back(temp_wake[4]);
				    } 
				  else{
				    throw
				      std::runtime_error("invalid specification of the wake type for 5 columns wake file");
				  } 
				  continue;
			    }
			    
			    if (num_columns==6) {
				  if (temp_wake.size()!=6)  throw std::runtime_error(" temp_wake size should be 6 in this case");				
				  z_coord.push_back(temp_wake[0]);				 				 
				  
				  if (get_wake_type()=="XLXTYLYTZ") {
				      xw_lead.push_back(temp_wake[1]);
				      xw_trail.push_back(temp_wake[2]);
				      yw_lead.push_back(temp_wake[3]);
				      yw_trail.push_back(temp_wake[4]);
				      z_wake.push_back(temp_wake[5]);
				  } 				
				  else{
				    throw
				      std::runtime_error("invalid specification of the wake type for 6 columns wake file");
				  } 
				  continue;
			    }
			    
			    throw
			    std::runtime_error("invalid specification of the wake type, the number of columnsin the wake file is too large "); 
			} 
		    } 
		}  // !line.empty()	    	    	    
	    }// while rfile
	    rfile.close();	  
	    size_wake=z_coord.size();
	    if (xw_lead.size()==0)   xw_lead.resize(size_wake, 0.0);
	    if (xw_trail.size()==0)  xw_trail.resize(size_wake, 0.0);
	    if (yw_lead.size()==0)   yw_lead.resize(size_wake, 0.0);
	    if (yw_trail.size()==0)  yw_trail.resize(size_wake, 0.0);
	    if ( z_wake.size()==0)   z_wake.resize(size_wake, 0.0);
	    
    
	    std::cout<<"  Wake_field: wake read from  "<<wake_file<<std::endl;
	    std::cout<<"  Wake_field: wake_type:  "<<get_wake_type()<<std::endl;
	  // wakes read!	 
    }//rank=0	  

  // Broadcasting to all     
	  int error=MPI_Bcast( (void *) &size_wake, 1, MPI_INT, 0,  MPI_COMM_WORLD );
	  z_coord.resize(size_wake);
	  xw_lead.resize(size_wake);
	  xw_trail.resize(size_wake);
	  yw_lead.resize(size_wake);
	  yw_trail.resize(size_wake);
	  z_wake.resize(size_wake);
	  
	  error=MPI_Bcast( z_coord.data(),  size_wake, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
	  error=MPI_Bcast( xw_lead.data(),  size_wake, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
	  error=MPI_Bcast( xw_trail.data(), size_wake, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
	  error=MPI_Bcast( yw_lead.data(),  size_wake, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
	  error=MPI_Bcast( yw_trail.data(), size_wake, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
	  error=MPI_Bcast( z_wake.data(),   size_wake, MPI_DOUBLE, 0,  MPI_COMM_WORLD );
	
	  if (z_coord[0]>0) {
	    double dz1=z_coord[1]-z_coord[0];
	    double dz2=z_coord[2]-z_coord[0];
	    delta_z=0.5*dz2-dz1;
	    istart=static_cast<int>(0.5*(1.-dz1/delta_z));
	    delta_z=0.25*dz2/(1.-istart);// delta_z recalcualted to reduce some numerical roundoff errors 
	    zstart=z_coord[0]-istart*istart*delta_z;
	  }
	  else if (z_coord[2]<=0) {
	    double dz1=z_coord[1]-z_coord[0];
	    double dz2=z_coord[2]-z_coord[0];
	    delta_z=dz1-0.5*dz2;
	    istart=static_cast<int>(0.5*(1.+dz1/delta_z));
	    delta_z=-0.25*dz2/(1.-istart);// delta_z recalcualted to reduce some numerical roundoff errors 
	    zstart=z_coord[0]+istart*istart*delta_z;
	  } 
	  else throw
	    std::runtime_error("wake file wrong: either the first z coordinate is positive or the first three  z  coordinates are negative");	  
  }
  catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;   
        MPI_Abort(MPI_COMM_WORLD, 454);
   }

}  

std::string Wake_field::get_wake_type() const { return wake_type;}
std::string Wake_field::get_wake_file_name() const { return wake_file;}
// std::vector<double> Wake_field::get_z_coord() const { return z_coord;}
// std::vector<double> Wake_field::get_xw_lead() const { return xw_lead;}
// std::vector<double> Wake_field::get_xw_trail() const { return xw_trail;}
// std::vector<double> Wake_field::get_yw_lead() const { return yw_lead;}
// std::vector<double> Wake_field::get_yw_trail() const{ return yw_trail;}
// std::vector<double> Wake_field::get_z_wake() const { return z_wake;}
int Wake_field::get_istart() const { return istart;}
double Wake_field::get_zstart() const { return zstart;}
double Wake_field::get_delta_z() const { return delta_z;}

MArray1d_ref const Wake_field::get_z_coord()  { 
  MArray1d_ref retval(z_coord.data(), boost::extents[z_coord.size()]);
  return retval;
}
MArray1d_ref const Wake_field::get_xw_lead() { 
  MArray1d_ref retval(xw_lead.data(), boost::extents[xw_lead.size()]);
  return retval;
} 

MArray1d_ref const Wake_field::get_xw_trail()  { 
  MArray1d_ref retval(xw_trail.data(), boost::extents[xw_trail.size()]);
  return retval;
} 

MArray1d_ref const Wake_field::get_yw_lead()  { 
  MArray1d_ref retval(yw_lead.data(), boost::extents[yw_lead.size()]);
  return retval;
}   

MArray1d_ref const Wake_field::get_yw_trail() { 
  MArray1d_ref retval(yw_trail.data(), boost::extents[yw_trail.size()]);
  return retval;
}   

MArray1d_ref const Wake_field::get_z_wake()  { 
  MArray1d_ref retval(z_wake.data(), boost::extents[z_wake.size()]);
  return retval;
}   




template<class Archive>
    void
    Wake_field::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(wake_file);
        ar & BOOST_SERIALIZATION_NVP(wake_type);
        ar & BOOST_SERIALIZATION_NVP(istart);
        ar & BOOST_SERIALIZATION_NVP(zstart);
	ar & BOOST_SERIALIZATION_NVP(delta_z);
	ar & BOOST_SERIALIZATION_NVP(z_coord);
	ar & BOOST_SERIALIZATION_NVP(xw_lead);
	ar & BOOST_SERIALIZATION_NVP(xw_trail);
	ar & BOOST_SERIALIZATION_NVP(yw_lead);
	ar & BOOST_SERIALIZATION_NVP(yw_trail);
	ar & BOOST_SERIALIZATION_NVP(z_wake);
    }


template
void
Wake_field::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Wake_field::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Wake_field::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Wake_field::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

