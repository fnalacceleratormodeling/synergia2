#include "fast_normal_form.h"
#include "Eigen/Eigen"
#include "synergia/foundation/physical_constants.h"

Fast_normal_form::Fast_normal_form(normalFormSage & nf):
  order(nf.order_),
  mappings_f(),
  mappings_g(),
  convert_from_cannon_sptr(),
  convert_to_canon_sptr(),
  closed_orbit(boost::extents[6]),
  eigenv(boost::extents[6][6]),
  inv_eigenv(boost::extents[6][6])
{
     for( int i = 0; i < order - 1; i++ ) {
            /// length=1. in Fast mapping constructors shoeld be irrelevant  
        mappings_f.push_back( CFast_mapping_sptr(new  CFast_mapping(nf.f_[i],1.)));    
        mappings_g.push_back( CFast_mapping_sptr(new  CFast_mapping(nf.g_[i],1.)));    
     }

     
   
    
    convert_from_cannon_sptr=Fast_mapping_sptr(new  Fast_mapping(nf.CanonToChef,1.));
    convert_to_canon_sptr=Fast_mapping_sptr(new  Fast_mapping(nf.ChefToCanon,1.));
    for (int i=0;i<6;++i){
      int chf_index=get_chef_index(i);
      closed_orbit[i]=nf.closed_orbit_(chf_index);
       for (int j=0;j<6;++j){
          eigenv[i][j]=nf.E_(i,j);
          inv_eigenv[i][j]=nf.invE_(i,j);
       }
    }

}


 CFast_mapping_sptr
 Fast_normal_form::get_map_f_sptr(int index) const
 {
   return mappings_f[index];
     
 }
 
 CFast_mapping_sptr
 Fast_normal_form::get_map_g_sptr(int index) const 
 {
   return mappings_g[index];
     
 }
 
 void 
 Fast_normal_form::convert_to_normal_form(Const_MArray1d_ref xyz,  MArray1dc_ref normal)
 {
   
   MArray1d coords(boost::extents[6]);
   for (int i=0;i<6;++i){       
       coords[i]=xyz[i]-closed_orbit[i];              
   }
 
   convert_to_canon_sptr->apply(coords);

   
   Eigen::VectorXcd u(6);
   Eigen::MatrixXcd ieigen(6,6);
   for (int i=0;i<6;++i){ 
      int chf_index=get_chef_index(i);
      u(chf_index)= coords[i];  
      for (int j=0;j<6;++j){
         ieigen(i,j)=inv_eigenv[i][j];
      }
   }
   
   u = ieigen*u;

   MArray1dc a(boost::extents[6]);
   std::complex<double> jj(0.,1.);

   for( int j = 0; j < 6; j++ ) {
             int chef_index=get_chef_index(j);              
             a[j]=u(chef_index);             
    }  

   for (int i = 0; i < order - 1; i++ ) { 
//        std::cout <<"fast i="<<i<< "    before: "<< std::endl ;
//        for( int j = 0; j < 6; j++ ) {     
//          std::cout << a[j]<< std::endl;
//        }
//       std::cout << std::endl;      
       mappings_f[i]->apply(a); 
//        for( int j = 0; j < 3; j++ ) { 
//            double a_real=0.5*real(a[2*j]+a[2*j+1]);
//            double a_img=0.5*imag(a[2*j]-a[2*j+1]);
//            a[2*j]=a_real+jj*a_img;
//            a[2*j+1]=a_real-jj*a_img;
//            
//        }
//        std::cout <<"fast i="<<i<< "    after: "<< std::endl ;
//        for( int j = 0; j < 6; j++ ) {        
//          std::cout << a[j] << std::endl;
//        } 
//         std::cout << std::endl;  
//         std::cout << std::endl;   
//        
       
   }
      
      
    normal=a;
        
   
 }

 void 
 Fast_normal_form::convert_from_normal_form(Const_MArray1dc_ref normal, MArray1d_ref xyz)
 {
   
   MArray1dc a(boost::extents[6]);
   for (int i=0;i<6;++i){       
      a[i]=normal[i];
   }
   
   
   for( int i = order - 2; i >= 0; --i ) {
//             std::cout <<"gfast i="<<i<< "    before: "<< std::endl ;
//             for( int j = 0; j < 6; j++ ) {     
//             std::cout << a[j]<< std::endl;
//             }
//             std::cout << std::endl;   
       mappings_g[i]->apply(a); 
//            std::cout <<"gfast i="<<i<< "    after: "<< std::endl ;
//             for( int j = 0; j < 6; j++ ) {        
//                 std::cout << a[j] << std::endl;
//             } 
//             std::cout << std::endl;  
//             std::cout << std::endl;              
//        
   }

   Eigen::VectorXcd u(6);
   Eigen::MatrixXcd eigen(6,6);
   for (int i=0;i<6;++i){ 
      int chf_index=get_chef_index(i);
      u(chf_index)= a[i];  
      for (int j=0;j<6;++j){
          eigen(i,j)=eigenv[i][j];
      }
    }
    u = eigen*u;

    for (int i = 0; i < 6; ++i ){
        int chf_index=get_chef_index(i);
        xyz[i]=real(u(chf_index));
        
    }

    convert_from_cannon_sptr->apply(xyz);

    for (int i=0;i<6;++i){       
       xyz[i] += closed_orbit[i];              
    } 
  
     
 }

 
void 
Fast_normal_form::convert_normal_to_xyz(MArray2d_ref coords)
{
    
   
    const MArray2d::size_type *coords_shape = coords.shape();
    const MArray2d::index *coords_bases = coords.index_bases();

    if ((coords_shape[1] != 7) || (coords_bases[1] != 0)) {
        throw std::runtime_error(
                "Fast_normal_form ::convert_normal_to_xyz expected nx[0:7] array");
    }
    for (unsigned int i = coords_bases[0];
            i != coords_bases[0] + coords_shape[0]; ++i) {
         MArray1dc a(boost::extents[6]);
         MArray1d  w(boost::extents[6]);
       
 
        for (int j = 0; j < 3; ++j) {
            a[2*j] = std::complex<double >(coords[i][2 * j],
                    coords[i][2 * j + 1]);
            a[2*j + 1] = std::conj(a[2*j]);
        }
        convert_from_normal_form(a, w);

        
        for (int j = 0; j < 6; ++j) {
            coords[i][j] = w[j];
        }
    } 
    
}   
    
void
Fast_normal_form::convert_xyz_to_normal(MArray2d_ref coords)
{
    
     
    const MArray2d::size_type *coords_shape = coords.shape();
    const MArray2d::index *coords_bases = coords.index_bases();



    if ((coords_shape[1] != 7) || (coords_bases[1] != 0)) {
        throw std::runtime_error(
                "Fast_normal_form::convert_xyz_to_normal expected nx[0:7] array");
    }

    for (unsigned int i = coords_bases[0];
            i != coords_bases[0] + coords_shape[0]; ++i) {
 
        
         MArray1d  w(boost::extents[6]);
         MArray1dc a(boost::extents[6]);
       
 
        for (int j = 0; j < 6; ++j) {
            w[j] = coords[i][j];
        }

        convert_to_normal_form(w, a);
       
 
        for (int j = 0; j < 3; ++j) {
            coords[i][2 * j] = a[2*j].real();
            coords[i][2 * j + 1] = a[2*j].imag();
        }
    }  
    
    
}    
 
 
std::vector<double>
Fast_normal_form::get_stationary_actions(const double stdx, const double stdy, const double stdct)
{
  

   Eigen::MatrixXd bmom(3,3);
   for (int i=0;i<3;++i){
     for (int j=0;j<3;++j){
        bmom(i,j)=2.0 * std::real( eigenv[i][j]*std::conj(eigenv[i][j]) ) ;
     }    
   }
   
   Eigen::MatrixXd inv_bmom(bmom.inverse());
   Eigen::VectorXd moments(3);
   moments(0) = stdx*stdx; 
   moments(1) = stdy*stdy;
   double stdt=stdct/pconstants::c;
   moments(2) = stdt*stdt;
   Eigen::VectorXd mact(inv_bmom*moments);
   
 
    std::vector<double> mean_actions(3);
    for (int i=0; i<3; ++i) {
       mean_actions[i] = mact(i);
    }

    return mean_actions;
  
}
    
template<class Archive>
    void
    Fast_normal_form::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(order);
        ar & BOOST_SERIALIZATION_NVP( mappings_f);
        ar & BOOST_SERIALIZATION_NVP( mappings_g);
        ar & BOOST_SERIALIZATION_NVP(convert_from_cannon_sptr); 
        ar & BOOST_SERIALIZATION_NVP(convert_to_canon_sptr);
        ar & BOOST_SERIALIZATION_NVP(closed_orbit);
        ar & BOOST_SERIALIZATION_NVP(eigenv);
        ar & BOOST_SERIALIZATION_NVP(inv_eigenv);         
    }


template
void
Fast_normal_form::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Fast_normal_form::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

// template
// void
// Fast_normal_form::serialize<boost::archive::text_oarchive >(
//         boost::archive::text_oarchive & ar, const unsigned int version);


template
void
Fast_normal_form::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Fast_normal_form::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);


// template
// void
// Fast_normal_form::serialize<boost::archive::text_iarchive >(
//         boost::archive::text_iarchive & ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(Fast_normal_form);