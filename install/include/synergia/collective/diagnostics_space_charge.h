#ifndef SP_DIAGNOSTICS_H_
#define SP_DIAGNOSTICS_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/deposit.h"

class Diagnostics_space_charge_rectangular : public Diagnostics
{ 
public:
  static const char name[];
  static const double field_fractional_beamsize_for_linear_interpolation;
private:
   
    bool have_writers;
    double s_n;    
    Hdf5_serial_writer<double > * writer_s_n;
    double num_macroparticles;
    Hdf5_serial_writer<double > * writer_num_macroparticles;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double s;
    Hdf5_serial_writer<double > * writer_s;  
    double step_length;
    Hdf5_serial_writer<double > * writer_step_length;  
    double factor_tune_Ederiv; // proportionality factor such as inch_tune_x= betax*factor_tune_Ederiv*dE/dx, inch_tune_y= betay*factor_tune_Ederiv*dE/dy 
                            // allows the possibility to extract the field derivatives dE/dx and dE/dy
    Hdf5_serial_writer<double > * writer_factor_tune_Ederiv; 
    MArray1d mean;
    Hdf5_serial_writer<MArray1d_ref > * writer_mean; 
    MArray1d std;
    Hdf5_serial_writer<MArray1d_ref > * writer_std;
    MArray1d step_betas;
    Hdf5_serial_writer<MArray1d_ref > * writer_step_betas;
    MArray2d inc_tune_shift;
    Hdf5_serial_writer<MArray2d_ref > * writer_inc_tune_shift;

    
public: 
    Diagnostics_space_charge_rectangular(std::string const& filename, std::string const& local_dir = "", 
                                         int const n_regression_points=4);
    Diagnostics_space_charge_rectangular();
    
           
    virtual bool
    is_serial() const;

    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

    virtual void
    update( );
    
    virtual void
    update(Bunch & bunch, Rectangular_grid const& En, int component, double time_step, double step_beta);

    virtual void
    write();
   
    MArray1d_ref 
    get_std() const; 
    
    double 
    get_num_macroparticles() const;
    
    double 
    get_step_length() const;
    
    double 
    get_factor_tune_Ederiv() const;
    
    MArray1d_ref 
    get_mean() const; 
    
    MArray1d_ref 
    get_step_betas() const; 
      
    MArray2d_ref 
    get_inc_tune_shift() const; 
    
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

     virtual
    ~Diagnostics_space_charge_rectangular();
    
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_space_charge_rectangular)
typedef boost::shared_ptr<Diagnostics_space_charge_rectangular > Diagnostics_space_charge_rectangular_sptr;
typedef std::list<Diagnostics_space_charge_rectangular_sptr >  Diagnostics_space_charge_rectangulars;


class Diagnostics_space_charge_3d_hockney : public Diagnostics
{ 
public:
  static const char name[];
  static const double field_fractional_beamsize_for_linear_interpolation;
private:
    bool have_writers;
    double s_n;    
    Hdf5_serial_writer<double > * writer_s_n;
    double num_macroparticles;
    Hdf5_serial_writer<double > * writer_num_macroparticles;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double s;
    Hdf5_serial_writer<double > * writer_s;    
    double step_length;
    Hdf5_serial_writer<double > * writer_step_length;  
    double factor_tune_Ederiv; // proportionality factor such as inch_tune_x= betax*factor_tune_Ederiv*dE/dx, inch_tune_y= betay*factor_tune_Ederiv*dE/dy 
                            // allows the possibility to extract the field derivatives dE/dx and dE/dy
    Hdf5_serial_writer<double > * writer_factor_tune_Ederiv; 
    MArray1d mean;
    Hdf5_serial_writer<MArray1d_ref > * writer_mean; 
    MArray1d std;
    Hdf5_serial_writer<MArray1d_ref > * writer_std;
    MArray1d step_betas;
    Hdf5_serial_writer<MArray1d_ref > * writer_step_betas;
    MArray2d inc_tune_shift;
    Hdf5_serial_writer<MArray2d_ref > * writer_inc_tune_shift;

    
public: 
    Diagnostics_space_charge_3d_hockney(std::string const& filename, std::string const& local_dir = "", 
                                        int const n_regression_points=4);
    Diagnostics_space_charge_3d_hockney();
    
    
    virtual bool
    is_serial() const;

    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

    virtual void
    update( );
    
    virtual void
    update(Bunch & bunch, Rectangular_grid const& En, int component, double time_step, double step_beta);

    virtual void
    write();
    
     MArray1d_ref 
    get_std() const; 
    
    double 
    get_num_macroparticles() const;
    
    double 
    get_step_length() const;
    
    double 
    get_factor_tune_Ederiv() const;
    
    MArray1d_ref 
    get_mean() const; 
    
    MArray1d_ref 
    get_step_betas() const; 
      
    MArray2d_ref 
    get_inc_tune_shift() const; 
   
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

     virtual
    ~Diagnostics_space_charge_3d_hockney();
    
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_space_charge_3d_hockney)
typedef boost::shared_ptr<Diagnostics_space_charge_3d_hockney > Diagnostics_space_charge_3d_hockney_sptr;
typedef std::list<Diagnostics_space_charge_3d_hockney_sptr >  Diagnostics_space_charge_3d_hockneys;






// class Diagnostics_space_charge_2d_bassetti_erskine : public Diagnostics
// {
//   public:
//   static const char name[];
// private:
//     bool have_writers;
//     double s_n;    
//     Hdf5_serial_writer<double > * writer_s_n;
//     double num_macroparticles;
//     Hdf5_serial_writer<double > * writer_num_macroparticles;
//     int repetition;
//     Hdf5_serial_writer<int > * writer_repetition;
//     double s;
//     Hdf5_serial_writer<double > * writer_s;    
//    // MArray1d mean;
//    // Hdf5_serial_writer<MArray1d_ref > * writer_mean; 
//     MArray1d std;
//     Hdf5_serial_writer<MArray1d_ref > * writer_std;
//     MArray1d step_betas;
//     Hdf5_serial_writer<MArray1d_ref > * writer_step_betas;
//    // MArray2d inc_tune_shift;
//    // Hdf5_serial_writer<MArray2d_ref > * writer_inc_tune_shift;
// 
//   public: 
//     Diagnostics_space_charge_2d_bassetti_erskine(std::string const& filename, std::string const& local_dir = "");
//     Diagnostics_space_charge_2d_bassetti_erskine();
//     
//     
//     virtual bool
//     is_serial() const;
// 
//     virtual void
//     init_writers(Hdf5_file_sptr file_sptr);
// 
//     virtual void
//     update( );
//     
//     virtual void
//     update(Bunch & bunch, double sigma_x, double sigma_y, double sigma_z, double time_step, double step_betax, double step_betay);
// 
//     virtual void
//     write();
//    
//     template<class Archive>
//         void
//         serialize(Archive & ar, const unsigned int version);
// 
//      virtual
//     ~Diagnostics_space_charge_2d_bassetti_erskine();    
//   
// };
// BOOST_CLASS_EXPORT_KEY(Diagnostics_space_charge_2d_bassetti_erskine)
// typedef boost::shared_ptr<Diagnostics_space_charge_2d_bassetti_erskine > Diagnostics_space_charge_2d_bassetti_erskine_sptr;
// typedef std::list<Diagnostics_space_charge_2d_bassetti_erskine_sptr >  Diagnostics_space_charge_2d_bassetti_erskines;

 #endif /* SP_DIAGNOSTICS_H_ */
 
