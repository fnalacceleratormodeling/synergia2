#include "fftw_helper.h"

void
Fftw_helper::construct(int *shape_in, bool z_periodic)
{
    shape = Int3(shape_in);
    shape[0] = 2 * shape_in[0];
    shape[1] = 2 * shape_in[1];
    if (! z_periodic) {
        shape[2] = 2 * shape_in[2];
    }
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
    if ((lower_limit == 0) || (local_nx == 0)) {
        left_guard = 0;
    } else {
        left_guard = 1;
    }
    if ((upper_limit >= shape[0] / 2) || (local_nx == 0)){
        right_guard = 0;
    } else {
        right_guard = 1;
    }
    data = reinterpret_cast<fftw_real *>
           (malloc(max_local_size * sizeof(fftw_real)));
    workspace = reinterpret_cast<fftw_real *>
                (malloc(max_local_size * sizeof(fftw_real)));
}

Fftw_helper::Fftw_helper(Int3 shape_in, bool z_periodic)
{
    construct(shape_in, z_periodic);
}

Fftw_helper::Fftw_helper(std::vector<int> shape_in, bool z_periodic)
{
    construct(&shape_in[0], z_periodic);
}

int
Fftw_helper::lower()
{
    return lower_limit;
}

int
Fftw_helper::upper()
{
    return upper_limit;
}

int
Fftw_helper::guard_lower()
{
    return lower_limit - left_guard;
}

int
Fftw_helper::guard_upper()
{
    return upper_limit + right_guard;
}

int
Fftw_helper::offset()
{
    return left_guard;
}

size_t
Fftw_helper::local_size()
{
    return max_local_size;
}

Int3
Fftw_helper::padded_shape_real()
{
    return Int3(shape[0], shape[1], 2*(shape[2] / 2 + 1));
}

Int3
Fftw_helper::padded_shape_complex()
{
    return Int3(shape[0], shape[1], shape[2] / 2 + 1);
}

void
Fftw_helper::transform(Real_scalar_field &in, Complex_scalar_field &out)
{
    size_t complex_data_length = (upper() - lower()) * shape[1] * (shape[2] / 2 + 1);
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
Fftw_helper::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
    size_t complex_data_length = (upper() - lower()) * shape[1] * (shape[2] / 2 + 1);
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

Fftw_helper::~Fftw_helper()
{
    fftw_free(data);
    fftw_free(workspace);
}
