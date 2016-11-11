#include "diagnostics_space_charge.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include "interpolate_rectangular_zyx.h"
#include "interpolate_rectangular_xyz.h"

const char Diagnostics_space_charge_rectangular::name[] = "spch_rectangular_diagnostics";
const double Diagnostics_space_charge_rectangular::field_fractional_beamsize_for_linear_interpolation=0.5;

Diagnostics_space_charge_rectangular::Diagnostics_space_charge_rectangular(std::string const& filename, 
        std::string const& local_dir, int const n_regression_points) :
        Diagnostics_space_charge_rectangular::Diagnostics(Diagnostics_space_charge_rectangular::name, filename, local_dir), 
        have_writers(false), s_n(0), writer_s_n(0), num_macroparticles(0), writer_num_macroparticles(0), 
        repetition(0), writer_repetition(0), s(0), writer_s(0),
        step_length(0), writer_step_length(0),
        factor_tune_Ederiv(0), writer_factor_tune_Ederiv(0),
        mean(boost::extents[3]), writer_mean(0), std(boost::extents[3]), writer_std(0), 
        step_betas(boost::extents[2]), writer_step_betas(0),
        inc_tune_shift(boost::extents[2][n_regression_points]), writer_inc_tune_shift(0)       
{
} 

Diagnostics_space_charge_rectangular::Diagnostics_space_charge_rectangular() : 
have_writers(false)
{
}

bool
Diagnostics_space_charge_rectangular::is_serial() const
{
    return true;
}

double 
Diagnostics_space_charge_rectangular::get_num_macroparticles() const
{
  return num_macroparticles;
} 

double 
Diagnostics_space_charge_rectangular::get_step_length() const
{
  return step_length;
} 

double 
Diagnostics_space_charge_rectangular::get_factor_tune_Ederiv() const
{
  return factor_tune_Ederiv;
} 


MArray1d_ref 
Diagnostics_space_charge_rectangular::get_mean() const
{
  return mean;
}


MArray1d_ref 
Diagnostics_space_charge_rectangular::get_std() const
{
  return std;
}

MArray1d_ref 
Diagnostics_space_charge_rectangular::get_step_betas() const
{
  return step_betas;
}


MArray2d_ref 
Diagnostics_space_charge_rectangular::get_inc_tune_shift() const
{
  return inc_tune_shift;
}
    

void
Diagnostics_space_charge_rectangular::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
        writer_num_macroparticles= new Hdf5_serial_writer<double > (file_sptr, "num_macroparticles");
        writer_repetition = new Hdf5_serial_writer<int > (file_sptr,  "repetition");
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");    
        writer_step_length= new Hdf5_serial_writer<double > (file_sptr, "step_length");
        writer_factor_tune_Ederiv= new Hdf5_serial_writer<double > (file_sptr, "factor_tune_Ederiv");
        writer_mean = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "mean");
        writer_std = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "std");
        writer_step_betas = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "step_betas");
        writer_inc_tune_shift= new Hdf5_serial_writer<MArray2d_ref > (file_sptr, "inc_tune_shift");
        have_writers = true;
    }
}


void
Diagnostics_space_charge_rectangular::update()
{
}


void
Diagnostics_space_charge_rectangular::update(Bunch & bunch, Rectangular_grid const& En, int component, double delta_t, double step_beta)
{ 

    bunch.convert_to_state(Bunch::fixed_z_lab);

    s_n = bunch.get_reference_particle().get_s_n();
    num_macroparticles= bunch.get_total_num();
    repetition = bunch.get_reference_particle().get_repetition();
    s = bunch.get_reference_particle().get_s();   
    step_betas[component]=step_beta;

    MArray1d mean_temp(boost::extents[6]);
    mean_temp=Core_diagnostics::calculate_mean(bunch);
    mean[0] = mean_temp[0];
    mean[1] = mean_temp[2];
    mean[2] = mean_temp[4]; 
    MArray1d std_temp(boost::extents[6]);
    std_temp=Core_diagnostics::calculate_std(bunch,mean_temp);  
    std[0] = std_temp[0];
    std[1] = std_temp[2];
    std[2] = std_temp[4];

    
    
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double gamma=bunch.get_reference_particle().get_gamma();
    double beta=bunch.get_reference_particle().get_beta();
// unit_conversion: [kg m/s] to [Gev/c] 
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
// scaled p = p/p_ref
    double p_ref=bunch.get_reference_particle().get_momentum();

    factor_tune_Ederiv=unit_conversion * q * delta_t/
           (4.0*mconstants::pi*p_ref*gamma*gamma*beta);
    double factor = step_beta*factor_tune_Ederiv*En.get_normalization();
   
    step_length=delta_t*pconstants::c*beta;
   
   
    Rectangular_grid_domain & domain(*En.get_domain_sptr());
    MArray3d_ref grid_points(En.get_grid_points());          
    int nderivs=inc_tune_shift.shape()[1];

    double bunch_mean=mean[component];  
    double bunch_std=std[component];
    //double delta=domain.get_cell_size()[component];
    double delta=field_fractional_beamsize_for_linear_interpolation*bunch_std/nderivs;
    double delta_x2=0.;
    double  y_deltax=0.; 
    
   
    
    if (component==0) {                        
      for (int i=0;i<nderivs;++i){
        double delta_h=(i+1)*delta;
        delta_x2 += 2.*delta_h*delta_h;
        y_deltax +=interpolate_rectangular_xyz(bunch_mean+delta_h,mean[1],mean[2], domain,grid_points)*delta_h;              
        y_deltax +=-interpolate_rectangular_xyz(bunch_mean-delta_h,mean[1],mean[2], domain,grid_points)*delta_h;
        inc_tune_shift[component][i]=factor*y_deltax/delta_x2; 
        //std::cout<<"i="<<i<<" dE/dx="<<y_deltax/delta_x2<<" dE="<<y_deltax<<" dx="<<delta_x2<<" factor="<<factor_tune_Ederiv<< std::endl;
  //              double  delta_h=bunch_mean+(i-nderivs/2)*delta;
  //              inc_tune_shift[component][i]=interpolate_rectangular_xyz(delta_h, mean[1], mean[2], domain,grid_points);
      }           
    }
    else if (component==1) {
      for (int i=0;i<nderivs;++i){
        double delta_h=(i+1)*delta;
        delta_x2 += 2.*delta_h*delta_h;
        y_deltax +=interpolate_rectangular_xyz(mean[0],bunch_mean+delta_h,mean[2], domain,grid_points)*delta_h;              
        y_deltax +=-interpolate_rectangular_xyz(mean[0],bunch_mean-delta_h,mean[2], domain,grid_points)*delta_h;
        inc_tune_shift[component][i]=factor*y_deltax/delta_x2;  
  //              double delta_h=bunch_mean+(i-nderivs/2)*delta;
  //              inc_tune_shift[component][i]=interpolate_rectangular_xyz(mean[0],delta_h,  mean[2], domain,grid_points);
    //  std::cout<<"i="<<i<<" dE/dy="<<y_deltax/delta_x2<<" dE="<<y_deltax<<" dy="<<delta_x2<<" factor="<<factor_tune_Ederiv<< std::endl;
      }                          
    }                              
}  


void
Diagnostics_space_charge_rectangular::write()
{
    // All bunch quantities are written in the bunch fixed t frame
    if (get_bunch().get_comm().has_this_rank()){
    if (get_write_helper().write_locally()) {
        init_writers(get_write_helper().get_hdf5_file_sptr());
        writer_s_n->append(s_n);
        writer_num_macroparticles->append(num_macroparticles);
        writer_repetition->append(repetition);
        writer_s->append(s);      
        writer_step_length->append(step_length);
        writer_factor_tune_Ederiv->append(factor_tune_Ederiv);
        writer_mean->append(mean); 
        writer_std->append(std);
        writer_step_betas->append(step_betas);
        writer_inc_tune_shift->append(inc_tune_shift);
        get_write_helper().finish_write();
    }
    }
}

template<class Archive>
    void
    Diagnostics_space_charge_rectangular::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
        ar & BOOST_SERIALIZATION_NVP(have_writers);
        ar & BOOST_SERIALIZATION_NVP(s_n);        
        ar & BOOST_SERIALIZATION_NVP(writer_s_n); 
        ar & BOOST_SERIALIZATION_NVP(num_macroparticles);
        ar & BOOST_SERIALIZATION_NVP(writer_num_macroparticles);
        ar & BOOST_SERIALIZATION_NVP(repetition);
        ar & BOOST_SERIALIZATION_NVP(writer_repetition);
        ar & BOOST_SERIALIZATION_NVP(s);
        ar & BOOST_SERIALIZATION_NVP(writer_s); 
        ar & BOOST_SERIALIZATION_NVP(step_length);
        ar & BOOST_SERIALIZATION_NVP(writer_step_length);
        ar & BOOST_SERIALIZATION_NVP(factor_tune_Ederiv);
        ar & BOOST_SERIALIZATION_NVP(writer_factor_tune_Ederiv);
        ar & BOOST_SERIALIZATION_NVP(mean);
        ar & BOOST_SERIALIZATION_NVP(writer_mean);
        ar & BOOST_SERIALIZATION_NVP(std);
        ar & BOOST_SERIALIZATION_NVP(writer_std);       
        ar & BOOST_SERIALIZATION_NVP(step_betas);
        ar & BOOST_SERIALIZATION_NVP(writer_step_betas);
        ar & BOOST_SERIALIZATION_NVP(inc_tune_shift);
        ar & BOOST_SERIALIZATION_NVP(writer_inc_tune_shift); 
    }

template
void
Diagnostics_space_charge_rectangular::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_space_charge_rectangular::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_space_charge_rectangular::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_space_charge_rectangular::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_space_charge_rectangular::~Diagnostics_space_charge_rectangular()
{
    if (have_writers) {
        delete writer_inc_tune_shift;
        delete writer_step_betas;
        delete writer_std;
        delete writer_mean;
        delete writer_s;
        delete writer_step_length;
        delete writer_factor_tune_Ederiv;
        delete writer_repetition;
        delete writer_s_n;
        delete writer_num_macroparticles;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_space_charge_rectangular)






// const char Diagnostics_space_charge_3d_hockney::name[] = "spch_3d_hockney_diagnostics";
// const double Diagnostics_space_charge_3d_hockney::field_fractional_beamsize_for_linear_interpolation=0.5;
// 
// Diagnostics_space_charge_3d_hockney::Diagnostics_space_charge_3d_hockney(std::string const& filename, 
//         std::string const& local_dir, int const n_regression_points) :
//         Diagnostics_space_charge_3d_hockney::Diagnostics(Diagnostics_space_charge_3d_hockney::name, filename, local_dir), 
//         have_writers(false), writer_s_n(0), writer_num_macroparticles(0), writer_repetition(0), writer_s(0), 
//         mean(boost::extents[3]), writer_mean(0), std(boost::extents[3]), writer_std(0), 
//         step_betas(boost::extents[2]), writer_step_betas(0),
//         inc_tune_shift(boost::extents[2][n_regression_points]), writer_inc_tune_shift(0)       
// {
// } 
// 
// Diagnostics_space_charge_3d_hockney::Diagnostics_space_charge_3d_hockney() : 
// have_writers(false)
// {
// }
// 
// bool
// Diagnostics_space_charge_3d_hockney::is_serial() const
// {
//     return true;
// }
// 
// void
// Diagnostics_space_charge_3d_hockney::init_writers(Hdf5_file_sptr file_sptr)
// {
//     if (!have_writers) {
//         writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
//         writer_num_macroparticles= new Hdf5_serial_writer<double > (file_sptr, "num_macroparticles");
//         writer_repetition = new Hdf5_serial_writer<int > (file_sptr,  "repetition");
//         writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");              
//         writer_mean = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "mean");
//         writer_std = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "std");
//         writer_step_betas = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "step_betas");
//         writer_inc_tune_shift= new Hdf5_serial_writer<MArray2d_ref > (file_sptr, "inc_tune_shift");
//         have_writers = true;
//     }
// }
// 
// 
// void
// Diagnostics_space_charge_3d_hockney::update()
// {
// }
// 
// 
// void
// Diagnostics_space_charge_3d_hockney::update(Bunch & bunch, Rectangular_grid const& En,  int component, double delta_t, double step_beta)
// { 
//   // the fields and bunch std are measured in the fixed_z_lab frame  
//          if(bunch.get_state() != Bunch::fixed_z_lab) {
//             std::cout<<" bunch state is="<<bunch.get_state()<<std::endl;
//             throw std::runtime_error(" bunch should be in fixed_t_bunch frame here and is not!");
//          }
//          s_n = bunch.get_reference_particle().get_s_n();
//          num_macroparticles= bunch.get_total_num();
//          repetition = bunch.get_reference_particle().get_repetition();
//          s = bunch.get_reference_particle().get_s();   
//          MArray1d mean_temp(boost::extents[6]);
//          mean_temp=Core_diagnostics::calculate_mean(bunch);
//          mean[0] = mean_temp[Bunch::x];
//          mean[1] = mean_temp[Bunch::y];
//          mean[2] = mean_temp[Bunch::z]; 
//          MArray1d std_temp(boost::extents[6]);
//          std_temp=Core_diagnostics::calculate_std(bunch,mean_temp);  
//          std[0] = std_temp[Bunch::x];
//          std[1] = std_temp[Bunch::y];
//          std[2] = std_temp[Bunch::z];
//          
//         
//       
//          /// En is electric field in units of N/C
//         double q = bunch.get_particle_charge() * pconstants::e; // [C]
//         // delta_t_beam: [s] in beam frame
//        // double delta_t_beam = delta_t / bunch.get_reference_particle().get_gamma();
//        
//         double gamma=bunch.get_reference_particle().get_gamma();
//         double beta=bunch.get_reference_particle().get_beta();
//          // unit_conversion: [kg m/s] to [Gev/c]
//         double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
//         // scaled p = p/p_ref
//         double p_scale = 1.0 / bunch.get_reference_particle().get_momentum();
//         double factor = step_beta*unit_conversion * q * delta_t * En.get_normalization()
//             * p_scale/(4.0*mconstants::pi)/gamma; // note the 1/gamma factor in delta_t_beam
//   
//          Rectangular_grid_domain & domain(*En.get_domain_sptr());
//          MArray3d_ref grid_points(En.get_grid_points());          
//          int nderivs=inc_tune_shift.shape()[1];
//         
//          double bunch_mean=mean[component];  
//          double bunch_std=std[component];
//          step_betas[component]=step_beta*delta_t*pconstants::c*beta;//step_beta*step_length;
//          //double delta=domain.get_cell_size()[component];
//          double delta=field_fractional_beamsize_for_linear_interpolation*bunch_std/nderivs;
//          double delta_x2=0.;
//          double  y_deltax=0.; 
//                    
//          if (component==0) {                        
//            for (int i=0;i<nderivs;++i){
//               double delta_h=(i+1)*delta;
//               delta_x2 += 2.*delta_h*delta_h;
//               y_deltax +=interpolate_rectangular_zyx(bunch_mean+delta_h,mean[1],mean[2], domain,grid_points)*delta_h;              
//               y_deltax +=-interpolate_rectangular_zyx(bunch_mean-delta_h,mean[1],mean[2], domain,grid_points)*delta_h;
//               inc_tune_shift[component][i]=factor*y_deltax/delta_x2; 
// //              double  delta_h=bunch_mean+(i-nderivs/2)*delta; // for testing
// //              inc_tune_shift[component][i]=interpolate_3d_hockney_xyz(delta_h, mean[1], mean[2], domain,grid_points);// for testing
//            }           
//          }
//          else if (component==1) {
//            for (int i=0;i<nderivs;++i){
//               double delta_h=(i+1)*delta;
//               delta_x2 += 2.*delta_h*delta_h;
//              y_deltax +=interpolate_rectangular_zyx(mean[0],bunch_mean+delta_h,mean[2], domain,grid_points)*delta_h;              
//              y_deltax +=-interpolate_rectangular_zyx(mean[0],bunch_mean-delta_h,mean[2], domain,grid_points)*delta_h;
//              inc_tune_shift[component][i]=factor*y_deltax/delta_x2;  
// //              double delta_h=bunch_mean+(i-nderivs/2)*delta;
// //              inc_tune_shift[component][i]=interpolate_3d_hockney_xyz(mean[0],delta_h,  mean[2], domain,grid_points);
//            }                          
//          }                              
// }  
// 
// 
// void
// Diagnostics_space_charge_3d_hockney::write()
// {
//     // All bunch quantities are written in the bunch fixed t frame
//     if (get_bunch().get_comm().has_this_rank()){
//     if (get_write_helper().write_locally()) {
//         init_writers(get_write_helper().get_hdf5_file_sptr());
//         writer_s_n->append(s_n);
//         writer_num_macroparticles->append(num_macroparticles);
//         writer_repetition->append(repetition);
//         writer_s->append(s);       
//         writer_mean->append(mean); 
//         writer_std->append(std);
//         writer_step_betas->append(step_betas);
//         writer_inc_tune_shift->append(inc_tune_shift);
//         get_write_helper().finish_write();
//     }
//     }
// }
// 
// template<class Archive>
//     void
//     Diagnostics_space_charge_3d_hockney::serialize(Archive & ar, const unsigned int version)
//     {
//         ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
//         ar & BOOST_SERIALIZATION_NVP(have_writers);
//         ar & BOOST_SERIALIZATION_NVP(s_n);        
//         ar & BOOST_SERIALIZATION_NVP(writer_s_n); 
//         ar & BOOST_SERIALIZATION_NVP(num_macroparticles);
//          ar & BOOST_SERIALIZATION_NVP(writer_num_macroparticles);
//         ar & BOOST_SERIALIZATION_NVP(repetition);
//         ar & BOOST_SERIALIZATION_NVP(writer_repetition);
//         ar & BOOST_SERIALIZATION_NVP(s);
//         ar & BOOST_SERIALIZATION_NVP(writer_s);      
//         ar & BOOST_SERIALIZATION_NVP(mean);
//         ar & BOOST_SERIALIZATION_NVP(writer_mean);
//         ar & BOOST_SERIALIZATION_NVP(std);
//         ar & BOOST_SERIALIZATION_NVP(writer_std);       
//         ar & BOOST_SERIALIZATION_NVP(step_betas);
//         ar & BOOST_SERIALIZATION_NVP(writer_step_betas);
//         ar & BOOST_SERIALIZATION_NVP(inc_tune_shift);
//         ar & BOOST_SERIALIZATION_NVP(writer_inc_tune_shift); 
//     }
// 
// template
// void
// Diagnostics_space_charge_3d_hockney::serialize<boost::archive::binary_oarchive >(
//         boost::archive::binary_oarchive & ar, const unsigned int version);
// 
// template
// void
// Diagnostics_space_charge_3d_hockney::serialize<boost::archive::xml_oarchive >(
//         boost::archive::xml_oarchive & ar, const unsigned int version);
// 
// template
// void
// Diagnostics_space_charge_3d_hockney::serialize<boost::archive::binary_iarchive >(
//         boost::archive::binary_iarchive & ar, const unsigned int version);
// 
// template
// void
// Diagnostics_space_charge_3d_hockney::serialize<boost::archive::xml_iarchive >(
//         boost::archive::xml_iarchive & ar, const unsigned int version);
// 
// Diagnostics_space_charge_3d_hockney::~Diagnostics_space_charge_3d_hockney()
// {
//     if (have_writers) {
//         delete writer_inc_tune_shift;
//         delete writer_step_betas;
//         delete writer_std;
//         delete writer_mean;
//         delete writer_s;
//         delete writer_repetition;
//         delete writer_s_n;
//         delete writer_num_macroparticles;
//     }
// }
// BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_space_charge_3d_hockney)
// 
// 
// 
// 
// 
// 
// const char Diagnostics_space_charge_2d_bassetti_erskine::name[] = "spch_bassetti_erskine_diagnostics";
// 
// 
// Diagnostics_space_charge_2d_bassetti_erskine::Diagnostics_space_charge_2d_bassetti_erskine(std::string const& filename, 
//         std::string const& local_dir) :
//         Diagnostics_space_charge_2d_bassetti_erskine::Diagnostics(Diagnostics_space_charge_2d_bassetti_erskine::name, filename, local_dir), 
//         have_writers(false), writer_s_n(0), writer_num_macroparticles(0), writer_repetition(0), writer_s(0), 
//        //mean(boost::extents[3]), writer_mean(0), 
//         std(boost::extents[3]), writer_std(0), 
//         step_betas(boost::extents[2]), writer_step_betas(0)
//       //  inc_tune_shift(boost::extents[2][n_regression_points]), writer_inc_tune_shift(0)       
// {
// } 
// 
// Diagnostics_space_charge_2d_bassetti_erskine::Diagnostics_space_charge_2d_bassetti_erskine() : 
// have_writers(false)
// {
// }
// 
// bool
// Diagnostics_space_charge_2d_bassetti_erskine::is_serial() const
// {
//     return true;
// }
// 
// void
// Diagnostics_space_charge_2d_bassetti_erskine::init_writers(Hdf5_file_sptr file_sptr)
// {
//     if (!have_writers) {
//         writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
//         writer_num_macroparticles= new Hdf5_serial_writer<double > (file_sptr, "num_macroparticles");
//         writer_repetition = new Hdf5_serial_writer<int > (file_sptr,  "repetition");
//         writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");              
//       //  writer_mean = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "mean");
//         writer_std = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "std");
//         writer_step_betas = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "step_betas");
//      //   writer_inc_tune_shift= new Hdf5_serial_writer<MArray2d_ref > (file_sptr, "inc_tune_shift");
//         have_writers = true;
//     }
// }
// 
// 
// void
// Diagnostics_space_charge_2d_bassetti_erskine::update()
// {
// }
// 
// 
// void
// Diagnostics_space_charge_2d_bassetti_erskine::update(Bunch & bunch, 
//                 double sigma_x, double sigma_y, double sigma_z, double time_step, double step_betax, double step_betay)
// { 
//   
//         s_n = bunch.get_reference_particle().get_s_n();
//         num_macroparticles= bunch.get_total_num();
//         repetition = bunch.get_reference_particle().get_repetition();
//         s = bunch.get_reference_particle().get_s();   
//         std[0] =sigma_x;
//         std[1] =sigma_y;
//         std[2] =sigma_z;
//         step_betas[0]=step_betax;
//         step_betas[1]=step_betay;
//       
// }  
// 
// 
// void
// Diagnostics_space_charge_2d_bassetti_erskine::write()
// {
//     // All bunch quantities are written in the bunch fixed t frame
//     if (get_bunch().get_comm().has_this_rank()){
//     if (get_write_helper().write_locally()) {
//         init_writers(get_write_helper().get_hdf5_file_sptr());
//         writer_s_n->append(s_n);
//         writer_num_macroparticles->append(num_macroparticles);
//         writer_repetition->append(repetition);
//         writer_s->append(s);       
//        // writer_mean->append(mean); 
//         writer_std->append(std);
//         writer_step_betas->append(step_betas);
//       //  writer_inc_tune_shift->append(inc_tune_shift);
//         get_write_helper().finish_write();
//     }
//     }
// }
// 
// template<class Archive>
//     void
//     Diagnostics_space_charge_2d_bassetti_erskine::serialize(Archive & ar, const unsigned int version)
//     {
//         ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
//         ar & BOOST_SERIALIZATION_NVP(have_writers);
//         ar & BOOST_SERIALIZATION_NVP(s_n);        
//         ar & BOOST_SERIALIZATION_NVP(writer_s_n); 
//         ar & BOOST_SERIALIZATION_NVP(num_macroparticles);
//          ar & BOOST_SERIALIZATION_NVP(writer_num_macroparticles);
//         ar & BOOST_SERIALIZATION_NVP(repetition);
//         ar & BOOST_SERIALIZATION_NVP(writer_repetition);
//         ar & BOOST_SERIALIZATION_NVP(s);
//         ar & BOOST_SERIALIZATION_NVP(writer_s);      
//        // ar & BOOST_SERIALIZATION_NVP(mean);
// //        ar & BOOST_SERIALIZATION_NVP(writer_mean);
//         ar & BOOST_SERIALIZATION_NVP(std);
//         ar & BOOST_SERIALIZATION_NVP(writer_std);       
//         ar & BOOST_SERIALIZATION_NVP(step_betas);
//         ar & BOOST_SERIALIZATION_NVP(writer_step_betas);
//       //  ar & BOOST_SERIALIZATION_NVP(inc_tune_shift);
//        // ar & BOOST_SERIALIZATION_NVP(writer_inc_tune_shift); 
//     }
// 
// template
// void
// Diagnostics_space_charge_2d_bassetti_erskine::serialize<boost::archive::binary_oarchive >(
//         boost::archive::binary_oarchive & ar, const unsigned int version);
// 
// template
// void
// Diagnostics_space_charge_2d_bassetti_erskine::serialize<boost::archive::xml_oarchive >(
//         boost::archive::xml_oarchive & ar, const unsigned int version);
// 
// template
// void
// Diagnostics_space_charge_2d_bassetti_erskine::serialize<boost::archive::binary_iarchive >(
//         boost::archive::binary_iarchive & ar, const unsigned int version);
// 
// template
// void
// Diagnostics_space_charge_2d_bassetti_erskine::serialize<boost::archive::xml_iarchive >(
//         boost::archive::xml_iarchive & ar, const unsigned int version);
// 
// Diagnostics_space_charge_2d_bassetti_erskine::~Diagnostics_space_charge_2d_bassetti_erskine()
// {
//     if (have_writers) {
//         //delete writer_inc_tune_shift;
//         delete writer_step_betas;
//         delete writer_std;
//        // delete writer_mean;
//         delete writer_s;
//         delete writer_repetition;
//         delete writer_s_n;
//         delete writer_num_macroparticles;
//     }
// }
// BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_space_charge_2d_bassetti_erskine)
