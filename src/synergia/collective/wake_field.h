#ifndef WAKE_FIELD_H_
#define WAKE_FIELD _H_

#include <string>
#include <vector>
#include "mpi.h"
#include <boost/shared_ptr.hpp>

class Wake_field
{
private:
  std::string wake_file;
  std::string wake_type; 
  
  ///assume the  wake functions are stored using a quadratic grid
  /// z[i]= (i-istart)^2*delta_z+zstart for i>istart
  /// z[i]= -(i-istart)^2*delta_z+zstart  for i< istart
  int istart;
  double zstart;
  double delta_z;
    
    
   std::vector<double> z_coord;    
   std::vector<double> xw_lead;// i.e. wake term proportional with the displacement of the leading (source) particle
   std::vector<double> xw_trail; // i.e. wake term proportional with the displacement of the trailing (affected) particle
   std::vector<double> yw_lead;// i.e. wake term proportional with the displacement of the leading (source) particle
   std::vector<double> yw_trail;// i.e. wake term proportional with the displacement of the trail particle
   std::vector<double> z_wake;
  

public:
  Wake_field();
  
  Wake_field(std::string const & wake_file, std::string const & wake_type);
  
  std::string get_wake_type() const;
  std::string get_wake_file_name() const;
  std::vector<double> get_z_coord() const;
  std::vector<double> get_xw_lead() const;
  std::vector<double> get_xw_trail() const; 
  std::vector<double> get_yw_lead() const;
  std::vector<double> get_yw_trail() const; 
  std::vector<double> get_z_wake() const;
  
}; 




typedef boost::shared_ptr<Wake_field> Wake_field_sptr;

#endif /* WAKE_FIELD_H_ */
