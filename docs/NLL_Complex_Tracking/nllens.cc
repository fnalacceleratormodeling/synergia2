#include <stdexcept>
#include <cmath>
#include <complex>

std::complex<double> Fpotential(const double x, const double y);
std::complex<double> Fderivative(const double x, const double y);
std::complex<double> carcsin(const std::complex<double> z);
std::complex<double> croot(const std::complex<double> z);

void NonlinearLensPropagatorCmplx(const double& knll, const double& cnll, double *coord)
{
/*
 !*****************************************************************
    ! The following subroutine computes the nonlinear momentum kick
    ! across a thin lens associated with a single short segment of the
    ! nonlinear magnetic insert described in V. Danilov and S. Nagaitsev,
    ! PRSTAB 13, 084002 (2010), Sect. V.A.  The arguments are as follows:
    !         knll - integrated strength of the lens (m)
    !         cnll - distance of singularities from the origin (m)
    !         coord = (x [m], px/p0, y [m], py/p0)
    ! This implementation is based on expressions in "Nonlinear Lens 
    ! Tracking in the IOTA Complex Potential," C. Mitchell, Feb. 9, 2017.
    ! Variable definitions are chosen to be consistent with TRACK_EXT.
    ! C. Mitchell 2/9/2017
    !*****************************************************************
*/

    std::complex<double> dF;
    std::complex<double> zeta;
    double x, y, kick, dPx, dPy;
    
    x = coord[0]/cnll;
    y = coord[2]/cnll;
    kick = -knll/cnll;
    // avoid the branch cuts
    if ((y == 0.0) && (std::abs(x) >= 1.0)) {
        throw std::runtime_error("nonlinear potential hits invalid branch cut region");
    }
    dF = Fderivative(x,y);
    dPx = kick*dF.real();
    dPy = -kick*dF.imag();
    
    coord[1] += dPx;
    coord[3] += dPy;
}

void InvariantPotentials(const double x, const double y, double& Hinv, double& Iinv)
{
/*
   !*****************************************************************
   ! This subroutine computes the dimensionless potentials that 
   ! define the contribution of the NLI to the two invariants 
   ! (H,I) for the IOTA ring.
   ! The arguments are as follows:
   !       (x,y) - normalized dimensionless coordinates
   !       Hinv - the vector potential describing the spatial
   !              dependence of the first invariant H.
   !       Iinv - function describing the spatial dependence
   !              of the second invariant I.
   ! This implementation is based on expressions in "Nonlinear Lens 
   ! Tracking in the IOTA Complex Potential," C. Mitchell, Feb. 9, 2017.
   ! C. Mitchell, 2/9/2017.
   !*****************************************************************
*/
    
    std::complex<double> zeta(x,y);
    std::complex<double> zetaconj = std::conj(zeta);
    std::complex<double> Hpotential = zeta/croot(zeta);
    std::complex<double> Ipotential = (zeta+zetaconj)/croot(zeta);
    Hpotential = Hpotential*carcsin(zeta);
    Ipotential = Ipotential*carcsin(zeta);
    Hinv = std::real(Hpotential);
    Iinv = std::real(Ipotential);
}

std::complex<double> Fpotential(const double x, const double y)
{
/*
   !****************************************
   ! Computes the dimensionless complex
   ! potential Az+i*Psi of the IOTA
   ! nonlinear insert.
   !****************************************
*/
    std::complex<double> zeta(x, y);
    std::complex<double> result = zeta/croot(zeta);
    result = result*carcsin(zeta);
    return(result);
}

std::complex<double> Fderivative(const double x, const double y)
{
/*
   !****************************************
   ! Computes the derivative of the
   ! dimensionless complex potential for
   ! the IOTA nonlinear insert.
   !****************************************
*/
    std::complex<double> zeta(x,y);
    std::complex<double> denom = croot(zeta);
    std::complex<double> result = zeta/(denom*denom);
    result = result + carcsin(zeta)/(denom*denom*denom);
    return(result);
}

std::complex<double> carcsin(const std::complex<double> z)
{
/*
   !******************************************
   ! Computes the complex function arcsin(z)
   ! using the principal branch.
   !******************************************
*/
    std::complex<double> c_i(0.0, 1.0);
    std::complex<double> result = c_i*z + croot(z);
    result = -c_i*std::log(result);
    return(result);
}

std::complex<double> croot(const std::complex<double> z)
{
/*
   !*******************************************
   ! Computes the complex function sqrt(1-z^2)
   ! using the principal branch.
   !*******************************************
*/
    std::complex<double> c_1(1.0, 0.0);
    return std::sqrt( c_1 - z*z );
}
