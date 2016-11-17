#ifndef EXCITE_MODE_H_
#define EXCITE_MODE_H_

#include "synergia/simulation/propagate_actions.h"
#include "synergia/collective/rectangular_grid.h"

class Excite_mode_actions : public Propagate_actions
{
  private:
    int turn_number_for_action;  
    
  public:
    Excite_mode_actions();
    Excite_mode_actions(int turn_number_for_action);
    virtual void
    turn_end_action(Stepper & stepper, Bunch_train & bunch_train, int turn_num);
    

    virtual void
    excite_bunch(Bunch &bunch);
                                     
   template<class Archive>
   void
   serialize(Archive & ar, const unsigned int version);
    
    virtual
   ~Excite_mode_actions();
    
};

BOOST_CLASS_EXPORT_KEY(Excite_mode_actions);
typedef boost::shared_ptr<Excite_mode_actions > Excite_mode_actions_sptr;


class Longitudinal_mode_actions : public Excite_mode_actions
{  
  private:
    MArray2d one_turn_map;
    double alpha_z;
    double beta_z;
    double delta_rz;
    int l_number;    
  public:
    Longitudinal_mode_actions();
    Longitudinal_mode_actions(MArray2d &one_turn_map, int turn_number_for_action, int l_number, double delta_rz=0.1);

    virtual void
    excite_bunch(Bunch &bunch);
  
   template<class Archive>
   void
   serialize(Archive & ar, const unsigned int version);
    
    virtual
   ~Longitudinal_mode_actions();
};

BOOST_CLASS_EXPORT_KEY(Longitudinal_mode_actions);
typedef boost::shared_ptr<Longitudinal_mode_actions > Longitudinal_mode_actions_sptr;

class Transverse_mode_actions : public Excite_mode_actions
{  
  private:
    MArray2d one_turn_map;
    double alpha_z;
    double beta_z;
    double alpha_t;
    double beta_t;
    double delta_rz; // fraction of zrms for displacement 
    int l_number; 
    int n_radial;
    int t_number; //t_number=1 is dipole
    double delta_rt; //fraction of xrms(yrms) for displacement 
    bool x_transverse; // true for horizontal modes, false for vertical modes
    bool strong_spc; // true for strong space charge modes, n_radial is the mode number
    bool from_file; // true if the excited moded is read from an text file
    std::string  excitation_mode_file_name;
    Rectangular_grid_sptr rectangular_grid_sptr; // 2d grid for excitation mode from file
    void read_mode();
  public:
    Transverse_mode_actions();
    Transverse_mode_actions(MArray2d &one_turn_map, int turn_number_for_action, bool x_transverse,
                            int l_number, int t_number, int n_radial=0, bool strong_spc=false, bool from_file=false,  double delta_rz=3.,
                               double delta_rt=0.3, std::string const & excitation_mode_file_name="");
                               
    Transverse_mode_actions(MArray2d &one_turn_map, int turn_number_for_action, bool x_transverse, double delta=1.);
    virtual void
    excite_bunch(Bunch &bunch);
  
//    template<class Archive>
//    void
//    serialize(Archive & ar, const unsigned int version);
//    
   
    template<class Archive>
       void
       save(Archive & ar, const unsigned int version) const;
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version);
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
    virtual
   ~Transverse_mode_actions();
};

BOOST_CLASS_EXPORT_KEY(Transverse_mode_actions);
typedef boost::shared_ptr<Transverse_mode_actions > Transverse_mode_actions_sptr;


inline double
interpolate_rectangular_xy(double x, double y, double left_x, double left_y, double cell_x, double cell_y,
         double shape_x, double shape_y,  MArray3d_ref const& a)
{

  int ix, iy;
  double offx, offy;
  double scaled_location;
  scaled_location = (x - left_x) / cell_x - 0.5;
  ix = fast_int_floor(scaled_location);
  offx = scaled_location - ix;
  
   scaled_location = (y - left_y) / cell_y - 0.5;
   iy = fast_int_floor(scaled_location);
   offy = scaled_location - iy;
  
  
    double val;
    if ( ix<0 ) {
      ix=0;
      offx=0.;
    }
      
    if ( ix >= (shape_x -1) ) {
      ix= shape_x -2;
      offx=1.;
    }
    
    
    if ( iy<0 ) {
      iy=0;
      offy=0.;
    }
    if ( iy >= (shape_y -1) ) {
      iy= shape_y -2;
      offy=1.;
    }
    
    val = ((1.0 - offx) * (1.0 - offy)*a[ix][iy][0]+offx * (1.0 - offy)*a[ix + 1][iy][0]+
            (1.0 - offx) * offy  * a[ix][iy + 1][0] + offx * offy * a[ix + 1][iy + 1][0]);
    
    
   
   /* 
   if ((ix >=0) && (ix < shape_x -1) && (iy >=0) &&
                    (iy < shape_y -1)){
      val = ((1.0 - offx) * (1.0 - offy)*a[ix][iy][0]+offx * (1.0 - offy)*a[ix + 1][iy][0]+
            (1.0 - offx) * offy  * a[ix][iy + 1][0] + offx * offy * a[ix + 1][iy + 1][0]);
    
    
    } 
    else if (ix<0) {
         
     }*/
    return val;    
}        

#endif /* EXCITE_MODE_H_ */
