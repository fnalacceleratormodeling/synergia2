#include "space_charge_2d_kv.h"

#include <complex>
#include <vector>

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/complex_error_function.h"

const int Space_charge_2d_kv::longitudinal_gaussian;
const int Space_charge_2d_kv::longitudinal_uniform;

const std::complex<double > complex_1(1.0, 0.0);
const std::complex<double > complex_0(0.0, 0.0);
const std::complex<double > complex_i(0.0, 1.0);

Space_charge_2d_kv::Space_charge_2d_kv() :
        Collective_operator("space charge")
        , sigma_x(0.0), sigma_y(0.0), sigma_cdt(0.0)
        , longitudinal_distribution(longitudinal_uniform)
{
}

Space_charge_2d_kv *
Space_charge_2d_kv::clone()
{
    return new Space_charge_2d_kv(*this);
}

// sets sigma_x, sigma_y for the distribution.  The sigmas are used
// to determine the extent of the uniform distribution.  sigma = total_range/sqrt(12)
void
Space_charge_2d_kv::set_sigma(double sigma_x, double sigma_y,
        double sigma_cdt)
{
    this->sigma_x = sigma_x;
    this->sigma_y = sigma_y;
    this->sigma_cdt = sigma_cdt;
    return;
}

int
Space_charge_2d_kv::get_longitudinal()
{
    return longitudinal_distribution;
}

void
Space_charge_2d_kv::set_longitudinal(int long_flag)
{
    if ((long_flag != longitudinal_gaussian) &&
        (long_flag != longitudinal_uniform)) {
        throw std::runtime_error("unknown longitudinal setting in Space_charge_2d_kv");
    }
    longitudinal_distribution = long_flag;
    return;
}

std::vector<double >
Space_charge_2d_kv::unit_efield(double arg_x, double arg_y)
{
    std::vector<double > retvec(2);
    unit_efield(arg_x, arg_y, retvec[0], retvec[1]);
    return retvec;
}

// returns Ex, Ey field for a unit kv charge
//
// from Miguel A. Furman, Compact Complex Expressions for the Electric
//   Field of 2-D Elliptical Charge Distributions: LBL-34682, CBP Note 014,
//   PEP-II/AP Note 34-93
//
void
Space_charge_2d_kv::unit_efield(double arg_x, double arg_y,
        double & E_x, double & E_y)
{
    // std of uniform distribution of length L = L/sqrt(12).
    // sqrt(12)*std is the full width which is twice the radius for
    // the KV distribution.
    double a = std::sqrt(3.0) * sigma_x;
    double b = std::sqrt(3.0) * sigma_y;

    bool inside = (arg_x/a)*(arg_x/a) + (arg_y/b)*(arg_y/b) < 1.0;
    bool xneg = arg_x<0.0;
    bool yneg = arg_y<0.0;
    if (inside) {
        E_x = 4 * arg_x/(a*(a+b));
        E_y = 4 * arg_y/(b*(a+b));
    } else {
        // for the exterior, we have to reflect to the first quadrant
        if (xneg) arg_x *= -1;
        if (yneg) arg_y *= -1;
        std::complex<double> zbar(arg_x, -arg_y);
        std::complex<double> E = 4.0/(zbar + std::sqrt(zbar*zbar - a*a + b*b));
        E_x = E.real();
        E_y = E.imag();
        if (xneg) E_x *= -1;
        if (yneg) E_y *= -1;
    }
    return;
}

void
Space_charge_2d_kv::apply(Bunch & bunch, double delta_t,
        Step & step, int verbosity, Logger & logger)
{
    bunch.convert_to_state(Bunch::fixed_z_lab);

    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d std(Core_diagnostics::calculate_std(bunch, mean));
    set_sigma(std[Bunch::x], std[Bunch::y], std[Bunch::z]);

    // dp/p kick =
    //
    //   N r_p * (1/gamma**2) *        delta-t   *         (1/(beta*gamma)) * unit_E_field
    //              E-B cancellation                               ^
    //                                                             |
    //                                                       1/p to get dp/p
    //

    double beta = bunch.get_reference_particle().get_beta();
    double gamma = bunch.get_reference_particle().get_gamma();

    double factor = pconstants::rp * pconstants::c * delta_t /
            (gamma*gamma*gamma*beta);
    // set longitudinal density depending on the longitudinal flag
    double total_q = bunch.get_real_num()*bunch.get_particle_charge();
    double line_charge_density;
    if (longitudinal_distribution == longitudinal_uniform) {
        line_charge_density = total_q / bunch.get_z_period_length() ;
    }

    double mean_x=mean[Bunch::x];
    double mean_y=mean[Bunch::y]; 
    double mean_z=mean[Bunch::z]; 
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x]-mean_x;
        double y = bunch.get_local_particles()[part][Bunch::y]-mean_y;
        double z = bunch.get_local_particles()[part][Bunch::z]-mean_z;

        if (longitudinal_distribution == longitudinal_gaussian) {
            line_charge_density = total_q
                                         * exp(-z * z / (2.0 * sigma_cdt * sigma_cdt))
                                         / (sqrt(2.0 * mconstants::pi) * sigma_cdt * beta);
        }
        double E_x, E_y;
        unit_efield(x, y, E_x, E_y);
        bunch.get_local_particles()[part][Bunch::xp] += E_x * factor * line_charge_density;
        bunch.get_local_particles()[part][Bunch::yp] += E_y * factor * line_charge_density;
    }
}

template<class Archive>
    void
    Space_charge_2d_kv::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(sigma_x);
        ar & BOOST_SERIALIZATION_NVP(sigma_y);
        ar & BOOST_SERIALIZATION_NVP(sigma_cdt);
        ar & BOOST_SERIALIZATION_NVP(longitudinal_distribution);
    }

template
void
Space_charge_2d_kv::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Space_charge_2d_kv::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Space_charge_2d_kv::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Space_charge_2d_kv::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Space_charge_2d_kv::~Space_charge_2d_kv()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Space_charge_2d_kv)
