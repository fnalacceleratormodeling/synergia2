#include <iostream>
#include <fftw.h>
#include <complex>

#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

void
doit()
{
    const int N(16);
    std::complex<double > in[N], out[N];
    for (int i = 0; i < N; ++i) {
        in[i] = 1.1 * i;
    }
    fftw_plan p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_one(p, (fftw_complex *) in, (fftw_complex *) out);
    std::cout << "in            out\n";
    for (int i = 0; i < N; ++i) {
        std::cout << in[i] << "    " << out[i] << std::endl;
    }
    fftw_destroy_plan(p);
}

using namespace boost::python;

BOOST_PYTHON_MODULE(trivial_fftw)
{
    def("doit", &doit);
}

