#include "space_charge_2d_bassetti_erskine.h"

#include <complex>
#include <vector>

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/complex_error_function.h"

const std::complex<double > complex_1(1.0, 0.0);
const std::complex<double > complex_0(0.0, 0.0);
const std::complex<double > complex_i(0.0, 1.0);

Space_charge_2d_bassetti_erskine::Space_charge_2d_bassetti_erskine() :
        Collective_operator("space charge")
{
}

Space_charge_2d_bassetti_erskine *
Space_charge_2d_bassetti_erskine::clone()
{
    return new Space_charge_2d_bassetti_erskine(*this);
}

void
Space_charge_2d_bassetti_erskine::set_sigma(double sigma_x, double sigma_y,
        double sigma_cdt)
{
    this->sigma_x = sigma_x;
    this->sigma_y = sigma_y;
    this->sigma_cdt = sigma_cdt;
    const double round_tolerance = 1.0e-6;
    is_round = (std::abs((sigma_x - sigma_y) / (sigma_x + sigma_y))
            < round_tolerance);
}

std::vector<double >
Space_charge_2d_bassetti_erskine::normalized_efield(double arg_x, double arg_y)
{
    std::vector<double > retvec(2);
    normalized_efield(arg_x, arg_y, retvec[0], retvec[1]);
    return retvec;
}

void
Space_charge_2d_bassetti_erskine::normalized_efield(double arg_x, double arg_y,
        double & E_x, double & E_y)
{
    enum
    {
        ur, ul, lr, ll
    } quadrant;

    double x = arg_x;
    double y = arg_y;

    // Asymptotic limit ...
    if ((sigma_x == 0.0) && (sigma_y == 0.0)) {
        double r_squared = x * x + y * y;
        if (r_squared < 1.0e-20) {
            throw std::runtime_error(
                    "Space_charge_2d_bassetti_erskine::normalized_efield: r is too small");
        }
        E_x = x / r_squared;
        E_y = y / r_squared;
    } else {
        // Round beam limit ...
        if (is_round) {
            double r_squared = x * x + y * y;
            double mean_sigma_squared = 2.0 * sigma_x * sigma_y;
            // Test for small r .....
            const double sigma_scale = 1.0e-6;
            if (r_squared > sigma_scale * mean_sigma_squared) {
                double vol_fact = (1.0 - exp(-r_squared / mean_sigma_squared));
                E_x = vol_fact * x / r_squared;
                E_y = vol_fact * y / r_squared;
            } else {
                E_x = x / mean_sigma_squared;
                E_y = y / mean_sigma_squared;
            }
        } else {
            // Elliptic beam ...
            if (arg_x >= 0.0) {
                if (arg_y >= 0.0) {
                    quadrant = ur;
                    x = arg_x;
                    y = arg_y;
                } else {
                    quadrant = lr;
                    x = arg_x;
                    y = -arg_y;
                }
            } else {
                if (arg_y >= 0.0) {
                    quadrant = ul;
                    x = -arg_x;
                    y = arg_y;
                } else {
                    quadrant = ll;
                    x = -arg_x;
                    y = -arg_y;
                }
            }
            // Check for normal processing ...
            bool normal = sigma_x > sigma_y;
            double tmp1;
            if (!normal) {
                tmp1 = sigma_x;
                sigma_x = sigma_y;
                sigma_y = tmp1;
                tmp1 = x;
                x = y;
                y = tmp1;
            }
            // The calculation ...
            double ds = sqrt(2.0 * (sigma_x * sigma_x - sigma_y * sigma_y));
            std::complex<double > z1 = (x + complex_i * y) / ds;
            double rsigma = sigma_y / sigma_x;
            std::complex<double > z2 = (x * rsigma + complex_i * y / rsigma)
                    / ds;
            double rx = x / sigma_x;
            double ry = y / sigma_y;
            double r2 = rx * rx + ry * ry;
            std::complex<double > z = -complex_i
                    * (wofz(z1) - exp(-r2 / 2.0) * wofz(z2)) / ds;

            if (normal) {
                if (quadrant == ur) {
                    E_x = real(z);
                    E_y = -imag(z);
                } else if (quadrant == ul) {
                    E_x = -real(z);
                    E_y = -imag(z);
                } else if (quadrant == lr) {
                    E_x = real(z);
                    E_y = imag(z);
                } else { // quadrant == ll
                    E_x = -real(z);
                    E_y = imag(z);
                }
            } else {
                if (quadrant == ur) {
                    E_x = -imag(z);
                    E_y = real(z);
                } else if (quadrant == ul) {
                    E_x = imag(z);
                    E_y = real(z);
                } else if (quadrant == lr) {
                    E_x = -imag(z);
                    E_y = -real(z);
                } else { // quadrant == ll
                    E_x = imag(z);
                    E_y = -real(z);
                }
            }
        }
    }
    return;
}

void
Space_charge_2d_bassetti_erskine::apply(Bunch & bunch, double delta_t,
        Step & step, int verbosity, Logger & logger)
{
    bunch.convert_to_state(Bunch::fixed_t);

    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d std(Core_diagnostics::calculate_std(bunch, mean));
    sigma_x = std[Bunch::x];
    sigma_y = std[Bunch::y];
    sigma_cdt = std[Bunch::z];

    // $\delta \vec{p} = \vec{F} \delta t = q \vec{E} \delta t$
    // point charge
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    // total charge
    double q_total = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    // delta_t_beam: [s] in beam frame
    double delta_t_beam = delta_t / bunch.get_reference_particle().get_gamma();
    // unit_conversion: [N] = [kg m/s^2] to [Gev/c]
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
    // scaled p = p/p_ref
    double p_scale = 1.0 / bunch.get_reference_particle().get_momentum();
    // conversion from normalized_efield to E_n/lambda in [(V/m)/(C/m)]
    double E_conversion = 1.0 / (2.0*sqrt(mconstants::pi)*pconstants::epsilon0);
    double factor = unit_conversion * q * delta_t_beam * p_scale * E_conversion;

    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x];
        double y = bunch.get_local_particles()[part][Bunch::y];
        double z = bunch.get_local_particles()[part][Bunch::z];
        // csp: This line charge density works only for the gaussian charge
        //      distribution.
        double line_charge_density = q_total
                * exp(-z * z / (2.0 * sigma_cdt * sigma_cdt))
                / (sqrt(2.0 * mconstants::pi) * sigma_cdt);
        double E_x, E_y;
        normalized_efield(x, y, E_x, E_y);
        bunch.get_local_particles()[part][Bunch::xp] += E_x * factor * line_charge_density;
        bunch.get_local_particles()[part][Bunch::yp] += E_y * factor * line_charge_density;
    }
}

template<class Archive>
    void
    Space_charge_2d_bassetti_erskine::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(sigma_x);
        ar & BOOST_SERIALIZATION_NVP(sigma_y);
        ar & BOOST_SERIALIZATION_NVP(sigma_cdt);
        ar & BOOST_SERIALIZATION_NVP(is_round);
    }

template
void
Space_charge_2d_bassetti_erskine::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Space_charge_2d_bassetti_erskine::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Space_charge_2d_bassetti_erskine::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Space_charge_2d_bassetti_erskine::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Space_charge_2d_bassetti_erskine::~Space_charge_2d_bassetti_erskine()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Space_charge_2d_bassetti_erskine)
