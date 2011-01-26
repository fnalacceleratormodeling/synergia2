#include "space_charge_2d_bassetti_erskine.h"

#include <complex>
#include <vector>

#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"

typedef std::complex<double > Complex;

const Complex complex_1(1.0, 0.0);
const Complex complex_0(0.0, 0.0);
const Complex complex_i(0.0, 1.0);

#include "basic_toolkit/MathConstants.h"

class BasErs_field
{
private:
    double sigma[2];
public:
    BasErs_field(double* = 0
    /* pointer to an array containing
     sigmax and sigmay [m] */
    );
    BasErs_field(const BasErs_field&);
    ~BasErs_field();

    char useRound; // By default = 1
    // If 1: then round beam approximation
    // used when horizontal and vertical
    // sigmas approximately equal.

    std::vector<double >
    NormalizedEField(double x, double y);
    /* returns the "normalized" electric field
     in the rest frame of the bunch, in inverse
     meters.  To get the field [V/m], this must
     be multiplied by Q/(2 pi epsilon_o), where
     Q is the line density of charge [C/m] (in
     rest frame). */

    void
    GetSigma(double*);
};

double const SIGMA_LIMIT = 64.0;
double const SIGMA_ROUND = 0.1;

BasErs_field::BasErs_field(double* sigin)
{
    int i;

    for (i = 0; i < 2; i++)
        sigma[i] = sigin[i];

    useRound = 1;
}

BasErs_field::~BasErs_field()
{
}

std::vector<double >
BasErs_field::NormalizedEField(double arg_x, double arg_y)
{
    std::vector<double > retvec(3);
    char normal;
    std::complex<double > z;
    double x, y;
    double sigmaX, sigmaY, ds, meanSigma;
    std::complex<double > arg1, arg2;
    double tmp1, r;
    std::complex<double > retarg1, retarg2;
    enum
    {
        ur, ul, lr, ll
    } quadrant;

    x = arg_x;
    y = arg_y;
    sigmaX = sigma[0];
    sigmaY = sigma[1];

    //std::cout << "input " <<sigmaX<<" "<<sigmaY<<" "<<x<<" "<<y<<std::endl;

    // Asymptotic limit ...
    if ((sigmaX == 0.0) && (sigmaY == 0.0)) {
        r = x * x + y * y;
        if (r < 1.0e-20) {
            std::cerr << "\n";
            std::cerr << "*** ERROR ***                                 \n";
            std::cerr << "*** ERROR *** BasErs::NormalizedEField        \n";
            std::cerr << "*** ERROR *** Asymptotic limit                \n";
            std::cerr << "*** ERROR *** r seems too small.              \n";
            std::cerr << "*** ERROR ***                                 \n";
            std::exit(1);
        }
        retvec[0] = x / r;
        retvec[1] = y / r;
        retvec[2] = 0.0;
        return retvec;
    }

    // Round beam limit ...
    if (useRound) {
        if ((fabs((sigmaX - sigmaY) / (sigmaX + sigmaY)) < SIGMA_ROUND)
                || ((pow(x / sigmaX, 2.0) + pow(y / sigmaY, 2.0)) > SIGMA_LIMIT)) {
            r = x * x + y * y;
            meanSigma = 2.0 * sigmaX * sigmaY;
            // Test for small r .....
            if (r > 1.0e-6 * meanSigma) {
                r = (1.0 - exp(-r / meanSigma)) / r;
                retvec[0] = x * r;
                retvec[1] = y * r;
                retvec[2] = 0.0;
                return retvec;
            } else {
                retvec[0] = x / meanSigma;
                retvec[1] = y / meanSigma;
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
    if (!(normal = (sigmaX > sigmaY))) {
        tmp1 = sigmaX;
        sigmaX = sigmaY;
        sigmaY = tmp1;
        tmp1 = x;
        x = y;
        y = tmp1;
    }

    // The calculation ...
    ds = sqrt(2.0 * (sigmaX * sigmaX - sigmaY * sigmaY));
    arg1 = x / ds + complex_i * y / ds;
    r = sigmaY / sigmaX;
    arg2 = ((x * r) / ds) + complex_i * ((y / r) / ds);

    retarg1 = w(arg1);
    retarg2 = w(arg2);

    // Normalization ...
    r = x / sigmaX;
    r = r * r;
    tmp1 = y / sigmaY;
    r += tmp1 * tmp1;

    z = retarg1;
    z -= retarg2 * exp(-r / 2.0);
    z *= -complex_i * MATH_SQRTPI / ds;

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

Space_charge_2d_bassetti_erskine::Space_charge_2d_bassetti_erskine() :
    Collective_operator("space charge")
{
}

void
Space_charge_2d_bassetti_erskine::apply(Bunch & bunch, double time_step,
        Step & step)
{

    double sigma[2];

    Diagnostics diagnostics(bunch);
    sigma[0] = diagnostics.get_std()[Bunch::x];
    sigma[1] = diagnostics.get_std()[Bunch::y];
    double sigma_cdt = diagnostics.get_std()[Bunch::z];

    BasErs_field bas_ers_field(sigma);

    // Alex:
    // In the lab frame  (Delta p) = q*E_eff* (Delta t) =factor*Efield*tau
    //
    //    what is factor=?
    //
    //1) the  arc length tau=beta*c* (Delta t), so (Delta t)= tau/(beta*c)
    //
    //2)   q=p/Brho=PH_CNV_brho_to_p
    //    because p unit is [GeV/c], the charge is measured in  q=c*10e-9
    //
    //    3) in the bunch frame
    //     E'= 1/(2*pi*eps0) *lambda'*Efield
    //     where Efield=normalized field, see BasErs_field.h
    //
    //E' --electric field in the bunch frame
    //    E=gamma* E' --electric field in the lab frame
    //
    //    E_eff=E-beta*B=E-Beta^2*E= E/gamma^2=E'/gamma
    //
    //4) charge density transformation:
    //lambda=lambda'*gamma ===>lambda'=lambda/gamma
    //
    //    5) lambda= current/v=current/beta*c
    //    or lambda=mbs.bunch_np*PH_MKS_e/mbs.z_length, but we did not defined mbs.z_length yet...
    //
    //6) Keep in mind the units used in the code p=p/p_ref_total.....

    double beta = bunch.get_reference_particle().get_beta();
    double gamma = bunch.get_reference_particle().get_gamma();
    double p_ref_total = bunch.get_reference_particle().get_momentum();
    double within_one_sigma_charge = bunch.get_real_num()
            * bunch.get_particle_charge() * erf(1.0 / std::sqrt(2.0))
            * pconstants::e;
    double current = within_one_sigma_charge / (2 * sigma_cdt / pconstants::c);

    double factor = PH_CNV_brho_to_p * current / (2.0 * mconstants::pi
            * pconstants::epsilon0 * beta * pconstants::c); //point 2) and 5) above
    factor = factor / (beta * pconstants::c); //         point 1) above
    factor = factor / (gamma * gamma); // point  3) and 4) above
    factor *= 1.0 / p_ref_total; // point 6)

    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x];
        double y = bunch.get_local_particles()[part][Bunch::y];
        std::vector<double > e_field(bas_ers_field.NormalizedEField(x, y));
        bunch.get_local_particles()[part][Bunch::xp] += e_field[0] * time_step
                * factor;
        bunch.get_local_particles()[part][Bunch::yp] += e_field[1] * time_step
                * factor;
    }
}

Space_charge_2d_bassetti_erskine::~Space_charge_2d_bassetti_erskine()
{
}
