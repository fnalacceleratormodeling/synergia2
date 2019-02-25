#ifndef FAST_NORMAL_H_
#define  FAST_NORMAL_H_

#include "fast_mapping.h"
#include "physics_toolkit/normalFormSage.h"
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>

/*
#include "synergia/foundation/reference_particle.h"
#include "synergia/bunch/bunch.h"
#include "synergia/lattice/chef_utils.h"

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#include "mxyzptlk/Mapping.h"
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

#include <vector>
#include <list>
#include <fstream>
#include <string>

#include "synergia/utils/serialization.h"*/

class Fast_normal_form 
{
  private:
     int order;
     std::vector<CFast_mapping_sptr >  mappings_f; //f_
     std::vector<CFast_mapping_sptr >  mappings_g; //g_
     Fast_mapping_sptr convert_from_cannon_sptr;     // CanonToChef,  but  synergia index order
     Fast_mapping_sptr convert_to_canon_sptr;  // ChefToCanon,  but  synergia index order
     MArray1d closed_orbit;  // closed_orbit_ , but  synergia index order
     MArray2dc eigenv; // E_;
     MArray2dc inv_eigenv; // invE_;
     
 
    
  public:
    Fast_normal_form(normalFormSage & nf);
    Fast_normal_form() = default;
    
    CFast_mapping_sptr
    get_map_f_sptr( int index) const;
    
    CFast_mapping_sptr
    get_map_g_sptr(int index) const;
    
    std::vector<double>
    get_stationary_actions(const double stdx, const double stdy, const double stdct);
    
    void 
    convert_to_normal_form(Const_MArray1d_ref xyz, MArray1dc_ref normal);
    
    void 
    convert_from_normal_form(Const_MArray1dc_ref normal, MArray1d_ref xyz);
    
    /// normal form below is stored in a real array as a[2j]=real(a) and a[2j+1]=Imag(a)
    void 
    convert_normal_to_xyz(MArray2d_ref coords);
     /// normal form below is stored in a real array as a[2j]=real(a) and a[2j+1]=Imag(a)
    void
    convert_xyz_to_normal(MArray2d_ref coords);
   
   template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};


BOOST_CLASS_EXPORT_KEY(Fast_normal_form);
typedef boost::shared_ptr<Fast_normal_form>  Fast_normal_form_sptr;

#endif /* FAST_NORMAL_H_ */
