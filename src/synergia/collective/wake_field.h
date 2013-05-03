#ifndef WAKE_FIELD_H_
#define WAKE_FIELD_H_

#include <string>
#include <vector>
#include "mpi.h"
#include <boost/shared_ptr.hpp>
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/serialization.h"

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
//   std::vector<double> get_z_coord() const;
//   std::vector<double> get_xw_lead() const;
//   std::vector<double> get_xw_trail() const; 
//   std::vector<double> get_yw_lead() const;
//   std::vector<double> get_yw_trail() const; 
//   std::vector<double> get_z_wake() const;
  MArray1d_ref const get_z_coord();
  MArray1d_ref const get_xw_lead();
  MArray1d_ref const get_xw_trail(); 
  MArray1d_ref const get_yw_lead();
  MArray1d_ref const get_yw_trail(); 
  MArray1d_ref const get_z_wake();

  int get_istart() const;
  double get_zstart() const;
  double get_delta_z() const;
  
  void multiply_xw_lead(double mltp);
  void multiply_xw_trail(double mltp);
  void multiply_yw_lead(double mltp);
  void multiply_yw_trail(double mltp);
  void multiply_z_wake(double mltp);
  
  template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);  
  
   ~Wake_field();
}; 



typedef boost::shared_ptr<Wake_field> Wake_field_sptr;

#endif /* WAKE_FIELD_H_ */
