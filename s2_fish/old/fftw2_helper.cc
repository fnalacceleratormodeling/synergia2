#include "fftw2_helper.h"
#include <algorithm>

Fftw2_helper_nompi::Fftw2_helper_nompi(Real_scalar_field &rho):
        Fftw_helper(rho)
{
    shape = rho.get_points().get_shape();
    shape.scale(2);
    timer("misc");
    plan = rfftwnd_create_plan(3, shape.c_array(),
                               FFTW_FORWARD, FFTW_MEASURE);
    inv_plan = rfftwnd_create_plan(3, shape.c_array(),
                                   FFTW_BACKWARD, FFTW_MEASURE);
    upper_limit = shape[0];
    max_local_size = shape[0] * shape[1] * shape[2];
}

int
Fftw2_helper_nompi::lower()
{
    return 0;
}

int
Fftw2_helper_nompi::upper()
{
    return upper_limit;
}

int
Fftw2_helper_nompi::offset()
{
    return 0;
}

size_t
Fftw2_helper_nompi::local_size()
{
    return max_local_size;
}

Int3
Fftw2_helper_nompi::padded_shape_real()
{
    return shape;
}

Int3
Fftw2_helper_nompi::padded_shape_complex()
{
    return shape;
}

void
Fftw2_helper_nompi::transform(Real_scalar_field &in, Complex_scalar_field &out)
{
    rfftwnd_one_real_to_complex(plan,
                                in.get_points().get_base_address(),
                                reinterpret_cast<fftw_complex *>
                                (out.get_points().get_base_address()));
}

void
Fftw2_helper_nompi::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
    rfftwnd_one_complex_to_real(inv_plan,
                                reinterpret_cast<fftw_complex *>
                                (in.get_points().get_base_address()),
                                out.get_points().get_base_address());
}

Fftw2_helper_nompi::~Fftw2_helper_nompi()
{}

Fftw2_helper_mpi::Fftw2_helper_mpi(Real_scalar_field &rho):
        Fftw_helper(rho)
{
    shape = rho.get_points().get_shape();
    shape.scale(2);
    timer("misc");
    plan = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, 3, shape.c_array(),
                                   FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    inv_plan = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, 3, shape.c_array(),
                                       FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    int local_nx, local_ny_after_transpose, local_y_start_after_transpose;
    rfftwnd_mpi_local_sizes(plan, &local_nx, &lower_limit,
                            &local_ny_after_transpose,
                            &local_y_start_after_transpose,
                            &max_local_size);
    upper_limit = lower_limit + local_nx;
    // padding for guard grids
    if (lower_limit == 0) {
        left_guard = 0;
    } else {
        left_guard = 1;
    }
    if (upper_limit >= rho.get_points().get_shape()[0]) {
        right_guard = 0;
    } else {
        right_guard = 1;
    }
    data = reinterpret_cast<fftw_real *>
           (malloc(max_local_size * sizeof(fftw_real)));
    workspace = reinterpret_cast<fftw_real *>
                (malloc(max_local_size * sizeof(fftw_real)));
}

int
Fftw2_helper_mpi::lower()
{
    return lower_limit;
}

int
Fftw2_helper_mpi::upper()
{
    return upper_limit;
}

int
Fftw2_helper_mpi::guard_lower()
{
    return lower_limit - left_guard;
}

int
Fftw2_helper_mpi::guard_upper()
{
    return upper_limit + right_guard;
}

int
Fftw2_helper_mpi::offset()
{
    return left_guard;
}

size_t
Fftw2_helper_mpi::local_size()
{
    return max_local_size;
}

Int3
Fftw2_helper_mpi::padded_shape_real()
{
    return Int3(shape[0], shape[1], 2*(shape[2] / 2 + 1));
}

Int3
Fftw2_helper_mpi::padded_shape_complex()
{
    return Int3(shape[0], shape[1], shape[2] / 2 + 1);
}

void
Fftw2_helper_mpi::transform(Real_scalar_field &in, Complex_scalar_field &out)
{
    size_t complex_data_length = (upper() - lower())*shape[1]*(shape[2] / 2 + 1);
    memcpy(reinterpret_cast<void*>(data),
           reinterpret_cast<void*>(in.get_points().get_offset_base_address(lower())),
           //~ in.get_points().get_length()*sizeof(double));
           complex_data_length*2*sizeof(double));
    rfftwnd_mpi(plan, 1,
                data,
                workspace,
                FFTW_NORMAL_ORDER);
    memcpy(reinterpret_cast<void*>(out.get_points().get_offset_base_address(lower())),
           reinterpret_cast<void*>(data),
           //~ out.get_points().get_length()*sizeof(std::complex<double>));
           complex_data_length*sizeof(std::complex<double>));
}

void
Fftw2_helper_mpi::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
    size_t complex_data_length = (upper() - lower())*shape[1]*(shape[2] / 2 + 1);
    memcpy(reinterpret_cast<void*>(data),
           reinterpret_cast<void*>(in.get_points().get_offset_base_address(lower())),
           complex_data_length*sizeof(std::complex<double>));
    rfftwnd_mpi(inv_plan, 1,
                data,
                workspace,
                FFTW_NORMAL_ORDER);
    memcpy(reinterpret_cast<void*>(out.get_points().get_offset_base_address(lower())),
           reinterpret_cast<void*>(data),
           complex_data_length*2*sizeof(double));
}

Fftw2_helper_mpi::~Fftw2_helper_mpi()
{
    fftw_free(data);
    fftw_free(workspace);
}
