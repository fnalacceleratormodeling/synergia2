#include "fftw_helper.h"


Fftw_helper::Fftw_helper(Real_scalar_field &rho)
{
}

int
Fftw_helper::lower()
{
  std::cout << "jfa: badbad lower\n";
  return 0;
}

int
Fftw_helper::upper()
{
  std::cout << "jfa: badbad upper\n";
  return 0;
}

int
Fftw_helper::guard_lower()
{
  std::cout << "jfa: badbad guard_lower\n";
  return 0;
}

int
Fftw_helper::guard_upper()
{
  std::cout << "jfa: badbad guard_upper\n";
  return 0;
}

int
Fftw_helper::offset()
{
  std::cout << "jfa: badbad offset\n";
  return 0;
}

size_t
Fftw_helper::local_size()
{
  std::cout << "jfa: badbad local_size\n";
  return 0;
}

Int3
Fftw_helper::padded_shape_real()
{
  std::cout << "jfa: badbad padded_shape_real\n";
  return Int3(1,1,1);
}

Int3
Fftw_helper::padded_shape_complex()
{
  std::cout << "jfa: badbad padded_shape_complex\n";
  return Int3(1,1,1);
}

void 
Fftw_helper::transform(Real_scalar_field &in, Complex_scalar_field &out)
{
  std::cout << "jfa: badbad transform\n";
}
 
void 
Fftw_helper::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
  std::cout << "jfa: badbad inv_transform\n";
}

Fftw_helper::~Fftw_helper()
{
}

