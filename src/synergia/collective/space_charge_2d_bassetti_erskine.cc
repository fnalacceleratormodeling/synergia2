#include "space_charge_2d_bassetti_erskine.h"

#include <complex>
#include <vector>

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/complex_error_function.h"

typedef std::complex<double > Complex;

const Complex complex_1(1.0, 0.0);
const Complex complex_0(0.0, 0.0);
const Complex complex_i(0.0, 1.0);


double const sigma_round = 0.1;

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
Space_charge_2d_bassetti_erskine::set_sigma(double sigma_x, double sigma_y, double sigma_cdt)
{
    this->sigma_x = sigma_x;
    this->sigma_y = sigma_y;
    this->sigma_cdt = sigma_cdt;
    const double round_tolerance = 1.0e-6;
    use_round = (std::abs(sigma_x/sigma_y -1) < round_tolerance);
}

std::vector<double >
Space_charge_2d_bassetti_erskine::normalized_efield(double arg_x, double arg_y)
{
    std::vector<double > retvec(3);
    double x, y;
    bool normal;
    std::complex<double > z;
    double ds, mean_sigma_squared;
    std::complex<double > arg1, arg2;
    double tmp1, r_squared;
    std::complex<double > retarg1, retarg2;
    enum
    {
        ur, ul, lr, ll
    } quadrant;

    x = arg_x;
    y = arg_y;

    // Asymptotic limit ...
    if ((sigma_x == 0.0) && (sigma_y == 0.0)) {
        r_squared = x * x + y * y;
        if (r_squared < 1.0e-20) {
            throw std::runtime_error("Asymptotic limit r seems too small in Space_charge_2d_bassetti_erskine::normalized_efield.");
        }
        retvec[0] = x / r_squared;
        retvec[1] = y / r_squared;
        retvec[2] = 0.0;
        return retvec;
    }

    // Round beam limit ...
    if (use_round) {
        if (fabs((sigma_x - sigma_y) / (sigma_x + sigma_y)) < sigma_round) {
            r_squared = x * x + y * y;
            mean_sigma_squared = 2.0 * sigma_x * sigma_y;
            // Test for small r .....
            if (r_squared > 1.0e-6 * mean_sigma_squared) {
             	double vol_fact = (1.0 - exp(-r_squared / mean_sigma_squared));
            	retvec[0] = vol_fact * x/r_squared;
            	retvec[1] = vol_fact * y/r_squared;
            	retvec[2] = 0.0;
            	return retvec;
            } else {
                retvec[0] = x / mean_sigma_squared;
                retvec[1] = y / mean_sigma_squared;
                retvec[2] = 0.0;
                return retvec;
            }
        }
    }

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
    normal = sigma_x > sigma_y;
    if (!normal) {
        tmp1 = sigma_x;
        sigma_x = sigma_y;
        sigma_y = tmp1;
        tmp1 = x;
        x = y;
        y = tmp1;
    }

    // The calculation ...
    ds = sqrt(2.0 * (sigma_x * sigma_x - sigma_y * sigma_y));
    arg1 = x / ds + complex_i * y / ds;
    double r = sigma_y / sigma_x;
    arg2 = ((x * r) / ds) + complex_i * ((y / r) / ds);

    retarg1 = wofz(arg1);
    retarg2 = wofz(arg2);

    // Normalization ...
    r = x / sigma_x;
    r = r * r;
    tmp1 = y / sigma_y;
    r += tmp1 * tmp1;

    z = retarg1;
    z -= retarg2 * exp(-r / 2.0);
    z *= -complex_i * sqrt(mconstants::pi) / ds;

    // And return ...
    retvec[2] = 0.0;
    if (normal) {
        if (quadrant == ur) {
            retvec[0] = real(z);
            retvec[1] = -imag(z);
            return retvec;
        }
        if (quadrant == ul) {
            retvec[0] = -real(z);
            retvec[1] = -imag(z);
            return retvec;
        }
        if (quadrant == lr) {
            retvec[0] = real(z);
            retvec[1] = imag(z);
            return retvec;
        }
        if (quadrant == ll) {
            retvec[0] = -real(z);
            retvec[1] = imag(z);
            return retvec;
        }
    } else {
        if (quadrant == ur) {
            retvec[0] = -imag(z);
            retvec[1] = real(z);
            return retvec;
        }
        if (quadrant == ul) {
            retvec[0] = imag(z);
            retvec[1] = real(z);
            return retvec;
        }
        if (quadrant == lr) {
            retvec[0] = -imag(z);
            retvec[1] = -real(z);
            return retvec;
        }
        if (quadrant == ll) {
            retvec[0] = imag(z);
            retvec[1] = -real(z);
            return retvec;
        }
        // ??? Just a guess; check this!
    }

    return retvec; // This line should never be reached.
}

void
Space_charge_2d_bassetti_erskine::apply(Bunch & bunch, double delta_t,
        Step & step, int verbosity, Logger & logger)
{
    // jfa: we should really convert to fixed_t state here and adjust
    //      factor accordingly.
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
    double factor = unit_conversion * q * delta_t_beam * p_scale
            / (2.0 * mconstants::pi * pconstants::epsilon0);

    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x];
        double y = bunch.get_local_particles()[part][Bunch::y];
        double z = bunch.get_local_particles()[part][Bunch::z];
        // csp: This line charge density works only for the gaussian charge
        //      distribution.
        double line_charge_density = q_total * exp(-z * z /(2.0 * sigma_cdt
                * sigma_cdt)) / (sqrt(2.0 * mconstants::pi) * sigma_cdt);
        double factor2 = line_charge_density * factor;
        std::vector<double > e_field(normalized_efield(x, y));
        bunch.get_local_particles()[part][Bunch::xp] += e_field[0] * factor2;
        bunch.get_local_particles()[part][Bunch::yp] += e_field[1] * factor2;
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
        ar & BOOST_SERIALIZATION_NVP(use_round);
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
