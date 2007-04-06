
#include "BasErs_field.h"

extern Complex w( Complex );

double const SIGMA_LIMIT = 64.0;
double const SIGMA_ROUND = 0.1;


BasErs_field::BasErs_field( double* sigin )
{
  int i;

  for( i = 0; i < 2; i++ ) sigma[i] = sigin[i];

  useRound = 1;
}

BasErs_field::~BasErs_field()
{
}

std::vector<double> BasErs_field::NormalizedEField( double arg_x, double arg_y )
{
  std::vector<double>  retvec(3);
  char    normal;
  Complex z;
  double  x, y;
  double  sigmaX, sigmaY, ds, meanSigma;
  Complex arg1, arg2;
  double  tmp1,r;
  Complex retarg1, retarg2;
  enum    { ur, ul, lr, ll } quadrant;

  x = arg_x;  
  y = arg_y;
  sigmaX = sigma[0];
  sigmaY = sigma[1];

  //std::cout << "input " <<sigmaX<<" "<<sigmaY<<" "<<x<<" "<<y<<std::endl;

  // Asymptotic limit ...
  if( ( sigmaX == 0.0 ) && ( sigmaY == 0.0 ) ) {
    r = x*x + y*y;
    if( r < 1.0e-20 ) {
      std::cerr << "\n";
      std::cerr << "*** ERROR ***                                 \n";
      std::cerr << "*** ERROR *** BasErs::NormalizedEField        \n";
      std::cerr << "*** ERROR *** Asymptotic limit                \n";
      std::cerr << "*** ERROR *** r seems too small.              \n";
      std::cerr << "*** ERROR ***                                 \n";
      exit(1);
    }
    retvec[0] = x/r;
    retvec[1] = y/r;
    retvec[2] = 0.0;
    return retvec;
  }

  // Round beam limit ...
  if( useRound ) {
    if( ( fabs( ( sigmaX - sigmaY ) / ( sigmaX + sigmaY ) ) < SIGMA_ROUND ) ||
        ( ( pow( x/sigmaX, 2.0 ) + pow( y/sigmaY, 2.0 ) ) > SIGMA_LIMIT )
      ) {
      r = x*x + y*y;
      meanSigma = 2.0*sigmaX*sigmaY;
      // Test for small r .....
      if( r > 1.0e-6*meanSigma ) {
  	r = ( 1.0 - exp(-r/ meanSigma ) ) / r;
  	retvec[0] = x*r;
  	retvec[1] = y*r;
  	retvec[2] = 0.0;
  	return retvec;
      }
      else {
  	retvec[0] = x/meanSigma;
  	retvec[1] = y/meanSigma;
  	retvec[2] = 0.0;
  	return retvec;
      }
    } 
  }


  // Elliptic beam ...
  if( arg_x >= 0.0 ) {
    if( arg_y >= 0.0 ) {
      quadrant = ur;
      x = arg_x;  
      y = arg_y;
    }
    else {
      quadrant = lr;
      x = arg_x;  
      y = - arg_y;
    }
  }
  else {
    if( arg_y >= 0.0 ) {
      quadrant = ul;
      x = - arg_x;  
      y = arg_y;
    }
    else {
      quadrant = ll;
      x = - arg_x;  
      y = - arg_y;
    }
  }

  // Check for normal processing ...
  if( !( normal = ( sigmaX > sigmaY ) ) ) {
   tmp1   = sigmaX;
   sigmaX = sigmaY;
   sigmaY = tmp1;
   tmp1   = x;
   x      = y;
   y      = tmp1;
  }

  // The calculation ...
  ds = sqrt(2.0*(sigmaX*sigmaX - sigmaY*sigmaY));
  arg1 = x/ds + complex_i*y/ds;  
  r = sigmaY/sigmaX;
  arg2 = ((x*r)/ds) + complex_i*((y/r)/ds);

  retarg1 = w( arg1 );
  retarg2 = w( arg2 );

  // Normalization ...
  r    = x/sigmaX;
  r    = r*r;
  tmp1 = y/sigmaY;
  r   += tmp1*tmp1;

  z    = retarg1;
  z   -= retarg2 * exp( - r/2.0 );
  z   *= - complex_i * MATH_SQRTPI / ds;

  // And return ...
  retvec[2] = 0.0;
  if( normal ) {
    if( quadrant == ur ) {
      retvec[0] =   real(z);
      retvec[1] = - imag(z);
      return retvec;
    }
    if( quadrant == ul ) {
      retvec[0] = - real(z);
      retvec[1] = - imag(z);
      return retvec;
    }
    if( quadrant == lr ) {
      retvec[0] =   real(z);
      retvec[1] =   imag(z);
      return retvec;
    }
    if( quadrant == ll ) {
      retvec[0] = - real(z);
      retvec[1] =   imag(z);
      return retvec;
    }
  }
  else {
    if( quadrant == ur ) {
      retvec[0] = - imag(z);
      retvec[1] =   real(z);
      return retvec;
    }
    if( quadrant == ul ) {
      retvec[0] =   imag(z);
      retvec[1] =   real(z);
      return retvec;
    }
    if( quadrant == lr ) {
      retvec[0] = - imag(z);
      retvec[1] = - real(z);
      return retvec;
    }
    if( quadrant == ll ) {
      retvec[0] =   imag(z);
      retvec[1] = - real(z);
      return retvec;
    }
    // ??? Just a guess; check this!
  }

  return retvec; // This line should never be reached.
}
