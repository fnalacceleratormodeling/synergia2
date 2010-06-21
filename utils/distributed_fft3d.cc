#include <iostream>
#include <cstring>

#include "distributed_fft3d.h"

Distributed_fft3d::Distributed_fft3d(std::vector<int > const & shape,
        bool z_periodic, Commxx const& comm, int planner_flags,
        std::string const& wisdom_filename) :
    shape(shape)
{
    //n.b. : we aren't using the wisdom_filename yet
#ifdef USE_FFTW2
    plan = rfftwnd_mpi_create_plan(comm.get(), 3, shape.c_array(),
            FFTW_REAL_TO_COMPLEX, planner_flags);
    inv_plan = rfftwnd_mpi_create_plan(comm.get(), 3, shape.c_array(),
            FFTW_COMPLEX_TO_REAL, planner_flags);
    int local_nx, local_ny_after_transpose, local_y_start_after_transpose;
    rfftwnd_mpi_local_sizes(plan, &local_nx, &lower_limit,
            &local_ny_after_transpose,
            &local_y_start_after_transpose,
            &max_local_size);
    data = reinterpret_cast<fftw_real *>
    (malloc(max_local_size * sizeof(fftw_real)));
    workspace = reinterpret_cast<fftw_real *>
    (malloc(max_local_size * sizeof(fftw_real)));
    if (max_local_size == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
    }
#else
    std::vector<int > padded_shape(get_padded_shape_real());
    fftw_mpi_init();
    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_3d(shape[0], shape[1],
            shape[2], comm.get(), &local_nx, &local_x_start);
    max_local_size = local_nx * padded_shape[1] * padded_shape[2];
    if (local_nx == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
        if (fftw_local_size > max_local_size) {
            max_local_size = fftw_local_size;
            std::cout << "jfa: warning: max_local_size had to be modified\n";
        }
    }
    data = (double *) fftw_malloc(sizeof(double) * max_local_size);
    workspace = (fftw_complex *) fftw_malloc(sizeof(double) * max_local_size);
    plan = fftw_mpi_plan_dft_r2c_3d(shape[0], shape[1], shape[2], data,
            workspace, comm.get(), planner_flags);
    inv_plan = fftw_mpi_plan_dft_c2r_3d(shape[0], shape[1], shape[2],
            workspace, data, comm.get(), planner_flags);
    lower_limit = local_x_start;
#endif //USE_FFTW2
    upper_limit = lower_limit + local_nx;
    // padding for guard grids
    if ((lower_limit == 0) || (local_nx == 0)) {
        left_guard = 0;
    } else {
        left_guard = 1;
    }
    if ((upper_limit >= shape[0] / 2) || (local_nx == 0)) {
        right_guard = 0;
    } else {
        right_guard = 1;
    }
}

int
Distributed_fft3d::get_lower() const
{
    return lower_limit;
}

int
Distributed_fft3d::get_upper() const
{
    return upper_limit;
}

int
Distributed_fft3d::get_guard_lower() const
{
    return lower_limit - left_guard;
}

int
Distributed_fft3d::get_guard_upper() const
{
    return upper_limit + right_guard;
}

int
Distributed_fft3d::get_offset() const
{
    return left_guard;
}

size_t
Distributed_fft3d::get_local_size() const
{
    return max_local_size;
}

std::vector<int >
Distributed_fft3d::get_padded_shape_real() const
{
    std::vector<int > retval(shape);
    retval[2] = 2 * (shape[2] / 2 + 1);
    return retval;
}

std::vector<int >
Distributed_fft3d::get_padded_shape_complex() const
{
    std::vector<int > retval(shape);
    retval[2] = shape[2] / 2 + 1;
    return retval;
}

void
Distributed_fft3d::transform(MArray3d_ref & in, MArray3dc_ref & out)
{
#ifdef USE_FFTW2
    size_t complex_data_length = (upper() - lower()) * shape[1] * (shape[2] / 2
            + 1);
    memcpy(reinterpret_cast<void* > (data),
            reinterpret_cast<void* > (in.origin() + in.index_bases()[0]
                    * in.strides()[0]), complex_data_length * 2
            * sizeof(double));
    rfftwnd_mpi(plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    memcpy(reinterpret_cast<void* > (out.origin() + out.index_bases()[0]
                    * out.strides()[0]), reinterpret_cast<void* > (data),
            complex_data_length * sizeof(std::complex<double >));
#else
    if (have_local_data) {
        memcpy(reinterpret_cast<void* > (data),
                reinterpret_cast<void* > (in.origin() + in.index_bases()[0]
                        * in.strides()[0]), get_local_size() * sizeof(double));
    }
    fftw_execute(plan);
    if (have_local_data) {
        memcpy(reinterpret_cast<void* > (out.origin() + out.index_bases()[0]
                * out.strides()[0]), reinterpret_cast<void* > (workspace),
                get_local_size() * sizeof(double));
    }
#endif //USE_FFTW2
}

void
Distributed_fft3d::inv_transform(MArray3dc_ref & in, MArray3d_ref & out)
{
#ifdef USE_FFTW2
    size_t complex_data_length = (upper() - lower()) * shape[1] * (shape[2] / 2
            + 1);
    memcpy(reinterpret_cast<void* > (data),
            reinterpret_cast<void* > (in.origin() + in.index_bases()[0]
                    * in.strides()[0]), complex_data_length
            * sizeof(std::complex<double >));
    rfftwnd_mpi(inv_plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    memcpy(reinterpret_cast<void* > (out.origin() + out.index_bases()[0]
                    * out.strides()[0]), reinterpret_cast<void* > (data),
            complex_data_length * 2 * sizeof(double));
#else
    if (have_local_data) {
        memcpy(reinterpret_cast<void* > (workspace),
                reinterpret_cast<void* > (in.origin() + in.index_bases()[0]
                        * in.strides()[0]), get_local_size() * sizeof(double));
    }
    fftw_execute(inv_plan);
    if (have_local_data) {
        memcpy(reinterpret_cast<void* > (out.origin() + out.index_bases()[0]
                * out.strides()[0]), reinterpret_cast<void* > (data),
                get_local_size() * sizeof(double));
    }
#endif //USE_FFTW2
}

Distributed_fft3d::~Distributed_fft3d()
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
