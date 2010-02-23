#include <iostream>
#include <cstring>

#include "fftw_helper.h"

void
Fftw_helper::construct(int *shape_in, bool z_periodic)
{
    shape = Int3(shape_in);
    shape[0] = 2 * shape_in[0];
    shape[1] = 2 * shape_in[1];
    if (! z_periodic) {
        shape[2] = 2 * shape_in[2];
    } else {shape[2] -= 1;
       }

    timer("misc");
#ifdef USE_FFTW2
    plan = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, 3, shape.c_array(),
                                   FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    inv_plan = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, 3, shape.c_array(),
                                       FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    int local_nx, local_ny_after_transpose, local_y_start_after_transpose;
    rfftwnd_mpi_local_sizes(plan, &local_nx, &lower_limit,
                            &local_ny_after_transpose,
                            &local_y_start_after_transpose,
                            &max_local_size);
    data = reinterpret_cast<fftw_real *>
           (malloc(max_local_size * sizeof(fftw_real)));
    workspace = reinterpret_cast<fftw_real *>
                (malloc(max_local_size * sizeof(fftw_real)));
#else
    Int3 padded_shape(padded_shape_real());
    fftw_mpi_init();
    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_3d(shape[0], shape[1],
            shape[2], MPI_COMM_WORLD, &local_nx, &local_x_start);
    max_local_size = local_nx * padded_shape[1] * padded_shape[2];
    if (fftw_local_size > max_local_size) {
        max_local_size = fftw_local_size;
        std::cout << "jfa: warning: max_local_size had to be modified\n";
    }
    data = (double *) fftw_malloc(sizeof(double) * max_local_size);
    workspace = (fftw_complex *) fftw_malloc(sizeof(double) * max_local_size);
    plan = fftw_mpi_plan_dft_r2c_3d(shape[0], shape[1], shape[2], data,
            workspace, MPI_COMM_WORLD, FFTW_ESTIMATE);
    inv_plan = fftw_mpi_plan_dft_c2r_3d(shape[0], shape[1], shape[2],
            workspace, data, MPI_COMM_WORLD, FFTW_ESTIMATE);
    lower_limit = local_x_start;
#endif //USE_FFTW2
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
#ifdef USE_FFTW2
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
#else
    memcpy(reinterpret_cast<void* > (data),
            reinterpret_cast<void* > (in.get_points().get_offset_base_address(
                    lower())), local_size() * sizeof(double));
    fftw_execute(plan);
    memcpy(reinterpret_cast<void* > (out.get_points().get_offset_base_address(
            lower())), reinterpret_cast<void* > (workspace), local_size()
            * sizeof(double));
#endif //USE_FFTW2
}

void
Fftw_helper::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
#ifdef USE_FFTW2
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
#else
    memcpy(reinterpret_cast<void* > (workspace),
            reinterpret_cast<void* > (in.get_points().get_offset_base_address(
                    lower())), local_size() * sizeof(double));
    fftw_execute(inv_plan);
    memcpy(reinterpret_cast<void* > (out.get_points().get_offset_base_address(
            lower())), reinterpret_cast<void* > (data), local_size()
            * sizeof(double));
#endif //USE_FFTW2
}

Fftw_helper::~Fftw_helper()
{
#ifdef USE_FFTW2
    fftw_free(data);
    fftw_free(workspace);
#else
    fftw_destroy_plan(plan);
    fftw_destroy_plan(inv_plan);
    fftw_free(data);
    fftw_free(workspace);
#endif //USE_FFTW2
}
