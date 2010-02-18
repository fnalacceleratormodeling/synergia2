#include <iostream>
// undefine symbols that conflict between iostream and mpich2
#if defined(SEEK_CUR)
#undef SEEK_CUR
#undef SEEK_SET
#undef SEEK_END
#endif /* defined(SEEK_CUR) */

#include <fftw_mpi.h>
#include <complex>

#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

void
doit()
{
    const int NX = 8, NY = 3;
    fftwnd_mpi_plan plan;
    std::complex<double > *data;

    plan = fftw2d_mpi_create_plan(MPI_COMM_WORLD, NX, NY, FFTW_FORWARD,
            FFTW_ESTIMATE);

    int local_nx, local_x_start;
    int local_ny_after_transpose, local_y_start_after_transpose;
    int total_local_size;
    fftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
            &local_ny_after_transpose, &local_y_start_after_transpose,
            &total_local_size);

    data = (std::complex<double >*) malloc(sizeof(std::complex<double >)
            * total_local_size);

    for (int x = 0; x < local_nx; ++x) {
        for (int y = 0; y < NY; ++y) {
            data[(x * NY + y)] = (x + local_x_start) * 10.0 + y * 0.1;
        }
    }

    for (int x = 0; x < local_nx; ++x) {
        for (int y = 0; y < NY; ++y) {
            std::cout << "pre:  " << x + local_x_start << "," << y << ": "
                    << data[(x * NY + y)] << std::endl;
            ;
        }
    }

    double t0, t1;
    t0 = MPI_Wtime();
    fftwnd_mpi(plan, 1, (fftw_complex*) data, NULL, FFTW_NORMAL_ORDER);
    t1 = MPI_Wtime();
    std::cout << "fft took " << t1-t0 << " secs\n";

    for (int x = 0; x < local_nx; ++x) {
        for (int y = 0; y < NY; ++y) {
            std::cout << "post: " << x + local_x_start << "," << y << ": "
                    << data[(x * NY + y)] << std::endl;
        }
    }

    fftwnd_mpi_destroy_plan(plan);
}

using namespace boost::python;

BOOST_PYTHON_MODULE(trivial_fftw)
{
    def("doit", &doit);
}

