#include <iostream>
#include <cstring>
#include <stdexcept>

#include "distributed_fft3d.h"

Distributed_fft3d::Distributed_fft3d(std::vector<int > const & shape,
        Commxx const& comm, int planner_flags,
        std::string const& wisdom_filename) :
    shape(shape), comm(comm), uppers(0), lengths(0)
{
    if (comm.get_size() / 2 >= shape[0] / 2) {
        throw std::runtime_error(
                "Distributed_fft3d: (number of processors)/2 must be <= shape[0]/2");
    }
    //n.b. : we aren't using the wisdom_filename yet
#ifdef USE_FFTW2
    throw std::runtime_error("Distributed_fft3d: fftw2 version needs updating and testing");
    plan = rfftwnd_mpi_create_plan(comm.get(), 3, shape.c_array(),
            FFTW_REAL_TO_COMPLEX, planner_flags);
    inv_plan = rfftwnd_mpi_create_plan(comm.get(), 3, shape.c_array(),
            FFTW_COMPLEX_TO_REAL, planner_flags);
    int local_nx, local_ny_after_transpose, local_y_start_after_transpose;
    rfftwnd_mpi_local_sizes(plan, &local_nx, &lower,
            &local_ny_after_transpose,
            &local_y_start_after_transpose,
            &local_size);
    data = reinterpret_cast<fftw_real *>
    (malloc(local_size * sizeof(fftw_real)));
    workspace = reinterpret_cast<fftw_real *>
    (malloc(local_size * sizeof(fftw_real)));
    if (local_size == 0) {
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
    local_size_real = local_nx * padded_shape[1] * padded_shape[2];
    int padded_shape2_complex = padded_shape[2] / 2;
    local_size_complex = local_nx * padded_shape[1] * padded_shape2_complex;
    local_size_allocated = local_size_real;
    if (local_nx == 0) {
        have_local_data = false;
        // create dummy space large enough for one local slice
        local_size_allocated = 1 * padded_shape[1] * padded_shape[2];
    } else {
        have_local_data = true;
        if (fftw_local_size > local_size_allocated) {
            std::cout
                    << "jfa: note: local_size_allocated had to be modified from "
                    << local_size_allocated << " to " << fftw_local_size
                    << " on rank " << comm.get_rank() << std::endl;
            local_size_allocated = fftw_local_size;
        }
    }
    data = (double *) fftw_malloc(sizeof(double) * local_size_allocated);
    workspace = (fftw_complex *) fftw_malloc(sizeof(double)
            * local_size_allocated);
    plan = fftw_mpi_plan_dft_r2c_3d(shape[0], shape[1], shape[2], data,
            workspace, comm.get(), planner_flags);
    inv_plan = fftw_mpi_plan_dft_c2r_3d(shape[0], shape[1], shape[2],
            workspace, data, comm.get(), planner_flags);
    lower = local_x_start;
#endif //USE_FFTW2
    upper = lower + local_nx;
}

Commxx &
Distributed_fft3d::get_comm()
{
    return comm;
}

int
Distributed_fft3d::get_lower() const
{
    return lower;
}

int
Distributed_fft3d::get_upper() const
{
    return upper;
}

void
Distributed_fft3d::calculate_uppers_lengths()
{
    if (uppers.empty()) {
        int size = comm.get_size();
        uppers.resize(size);
        lengths.resize(size);
        if (size == 1) {
            uppers.at(0) = upper;
            lengths.at(0) = upper - lower;
        } else {
            MPI_Allgather((void*) (&upper), 1, MPI_INT, (void*) (&uppers[0]),
                    1, MPI_INT, comm.get());
        }
        lengths[0] = uppers[0]*shape[1]*shape[2];
        for (int i = 1; i < size; ++i) {
            lengths[i] = (uppers[i] - uppers[i - 1]) * shape[1] * shape[2];
            if (lengths[i] < 0) {
                lengths[i] = 0;
            }
        }
    }
}

std::vector<int > const&
Distributed_fft3d::get_uppers()
{
    calculate_uppers_lengths();
    return uppers;
}

std::vector<int > const&
Distributed_fft3d::get_lengths()
{
    calculate_uppers_lengths();
    return lengths;
}

std::vector<int > const&
Distributed_fft3d::get_shape() const
{
    return shape;
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
            reinterpret_cast<void* > (in.origin() + (in.index_bases()[0] + lower)
                    * in.strides()[0]), complex_data_length * 2
            * sizeof(double));
    rfftwnd_mpi(plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    memcpy(reinterpret_cast<void* > (out.origin() + (out.index_bases()[0] + lower)
                    * out.strides()[0]), reinterpret_cast<void* > (data),
            complex_data_length * sizeof(std::complex<double >));
#else
    if (have_local_data) {
        memcpy(
                reinterpret_cast<void* > (data),
                reinterpret_cast<void* > (in.origin() + lower * in.strides()[0]),
                local_size_real * sizeof(double));
    }
    fftw_execute(plan);
    if (have_local_data) {
        memcpy(reinterpret_cast<void* > (out.origin() + lower
                * out.strides()[0]), reinterpret_cast<void* > (workspace),
                local_size_complex * sizeof(std::complex<double >));
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
            reinterpret_cast<void* > (in.origin() + (in.index_bases()[0] + lower)
                    * in.strides()[0]), complex_data_length
            * sizeof(std::complex<double >));
    rfftwnd_mpi(inv_plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    memcpy(reinterpret_cast<void* > (out.origin() + (out.index_bases()[0] + lower)
                    * out.strides()[0]), reinterpret_cast<void* > (data),
            complex_data_length * 2 * sizeof(double));
#else
    if (have_local_data) {
        memcpy(
                reinterpret_cast<void* > (workspace),
                reinterpret_cast<void* > (in.origin() + lower * in.strides()[0]),
                local_size_complex * sizeof(std::complex<double >));
    }
    fftw_execute(inv_plan);
    if (have_local_data) {
        memcpy(reinterpret_cast<void* > (out.origin() + lower
                * out.strides()[0]), reinterpret_cast<void* > (data),
                local_size_real * sizeof(double));
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
