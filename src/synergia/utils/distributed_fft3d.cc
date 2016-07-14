#include <iostream>
#include <cstring>
#include <stdexcept>
#include "multi_array_offsets.h"
#include "distributed_fft3d.h"

#include <omp.h>

Distributed_fft3d::Distributed_fft3d(std::vector<int > const & shape,
        Commxx_sptr comm_sptr, int planner_flags,
        std::string const& wisdom_filename) :
        uppers(0), lengths(0), shape(shape), comm_sptr(comm_sptr)
{
    if (comm_sptr->get_size() / 2 > shape[0] / 2) {
        throw std::runtime_error(
                "Distributed_fft3d: (number of processors)/2 must be <= shape[0]/2");
    }

    std::vector<int > padded_shape(get_padded_shape_real());
    int padded_shape2_complex = padded_shape[2] / 2;

    //n.b. : we aren't using the wisdom_filename yet
#ifdef USE_FFTW2
    plan = rfftwnd_mpi_create_plan(comm_sptr->get(), 3, &shape[0],
            FFTW_REAL_TO_COMPLEX, planner_flags);
    inv_plan = rfftwnd_mpi_create_plan(comm_sptr->get(), 3, &shape[0],
            FFTW_COMPLEX_TO_REAL, planner_flags);
    // fftw2 mpi manual says to allocate data and workspace of
    // size total_local_size, not local_nx size.
    // http://www.fftw.org/fftw2_doc/fftw_4.html#SEC47
    int total_local_size;
    int local_nx, local_x_start;
    int local_ny_after_transpose, local_y_start_after_transpose;
    rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
            &local_ny_after_transpose,
            &local_y_start_after_transpose,
            &total_local_size);
    local_size_real = local_nx * padded_shape[1] * padded_shape[2];
    local_size_complex = local_nx * padded_shape[1] * padded_shape2_complex;
    data = (fftw_real *) malloc(total_local_size * sizeof(fftw_real));
    workspace = (fftw_real *) malloc(total_local_size * sizeof(fftw_real));
    if (local_size_real == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
    }
#else
    // quoted from 
    // http://www.fftw.org/doc/Combining-MPI-and-Threads.html#Combining-MPI-and-Threads
    //
    //   "we must call fftw_init_threads before fftw_mpi_init. 
    //    This is critical for technical reasons having to do 
    //    with how FFTW initializes its list of algorithms."

    fftw_init_threads();
    fftw_mpi_init();

    int num_threads;
    #pragma omp parallel shared(num_threads)
    { num_threads = omp_get_num_threads(); }

    fftw_plan_with_nthreads(num_threads);

    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_3d(shape[0], shape[1],
            shape[2], comm_sptr->get(), &local_nx, &local_x_start);
    local_size_real = local_nx * padded_shape[1] * padded_shape[2];
    local_size_complex = local_nx * padded_shape[1] * padded_shape2_complex;
    if (local_nx == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
    }
    // MEDIUM HACK. fftw often (always?) returns a local size that is
    // impossibly small. Adjust to the smallest possible size.
    if (local_size_real > fftw_local_size) {
        fftw_local_size = local_size_real;
    }
    // GREAT BIG HACK. fftw sometimes crashes when memory_fudge_factor = 1. The test
    // program test_distributed_fft3d_mpi fails when memory_fudge_factor = 1 and
    // numproc = 3. This is either a bug in fftw3, or a failing of the documentation.
    // The precise value "2" is just a guess.
    const int memory_fudge_factor = 2;
    data = (double *) fftw_malloc(sizeof(double) * fftw_local_size
            * memory_fudge_factor);
    workspace = (fftw_complex *) fftw_malloc(sizeof(double) * fftw_local_size
            * memory_fudge_factor);
    plan = fftw_mpi_plan_dft_r2c_3d(shape[0], shape[1], shape[2], data,
            workspace, comm_sptr->get(), planner_flags);
    inv_plan = fftw_mpi_plan_dft_c2r_3d(shape[0], shape[1], shape[2],
            workspace, data, comm_sptr->get(), planner_flags);
#endif //USE_FFTW2
    lower = local_x_start;
    upper = lower + local_nx;
}

Commxx_sptr
Distributed_fft3d::get_comm_sptr()
{
    return comm_sptr;
}

Commxx &
Distributed_fft3d::get_comm()
{
    return *comm_sptr;
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
        int size = comm_sptr->get_size();
        uppers.resize(size);
        lengths.resize(size);
        if (size == 1) {
            uppers[0] = upper;
        } else {
            MPI_Allgather((void*) (&upper), 1, MPI_INT, (void*) (&uppers[0]),
                    1, MPI_INT, comm_sptr->get());
        }
        lengths[0] = uppers[0] * shape[1] * shape[2];
        for (int i = 1; i < size; ++i) {
            if (uppers[i] == 0) {
                uppers[i] = uppers[i - 1];
            }
            lengths[i] = (uppers[i] - uppers[i - 1]) * shape[1] * shape[2];
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
    if (have_local_data) {
        if (in.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible first index offset in input array");
        }
        if (static_cast<int >(in.index_bases()[0] + in.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != get_padded_shape_real()[1]) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible second dimension of input array");
        }
        if (static_cast<int >(in.shape()[2]) != get_padded_shape_real()[2]) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible third dimension of input array");
        }
        if (static_cast<int >(out.index_bases()[0]) > lower) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible first index offset in output array");
        }
        if (static_cast<int >(out.index_bases()[0] + out.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible first dimension of output array");
        }
        if (static_cast<int >(out.shape()[1])
                != get_padded_shape_complex()[1]) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible second dimension of output array");
        }
        if (static_cast<int >(out.shape()[2])
                != get_padded_shape_complex()[2]) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible third dimension of output array");
        }
    }
#ifdef USE_FFTW2
    if (have_local_data) {
        memcpy((void*) data,(void*) multi_array_offset(in, lower, 0, 0),
                local_size_real * sizeof(double));
    }
    rfftwnd_mpi(plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    if (have_local_data) {
        memcpy((void* ) multi_array_offset(out, lower, 0, 0),
                (void*) (data),
                local_size_complex * sizeof(std::complex<double >));
    }

#else
    if (have_local_data) {
        memcpy((void*) data, (void*) multi_array_offset(in, lower, 0, 0),
                local_size_real * sizeof(double));
    }
    fftw_execute(plan);
    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0, 0),
                (void*) (workspace), local_size_complex * sizeof(std::complex<
                        double >));
    }
#endif //USE_FFTW2
}


void
Distributed_fft3d::transform(MArray3d_ref & in, fftw_complex * out)
{
    if (have_local_data) {
        if (in.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible first index offset in input array");
        }
        if (static_cast<int >(in.index_bases()[0] + in.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != get_padded_shape_real()[1]) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible second dimension of input array");
        }
        if (static_cast<int >(in.shape()[2]) != get_padded_shape_real()[2]) {
            throw std::runtime_error(
                    "Distributed_fft3d::transform found an incompatible third dimension of input array");
        }
    }
#ifdef USE_FFTW2
    if (have_local_data) {
        memcpy((void*) data,(void*) multi_array_offset(in, lower, 0, 0),
                local_size_real * sizeof(double));
    }
    rfftwnd_mpi(plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    if (have_local_data) {
        memcpy((void* ) out, (void*) (data),
                local_size_complex * sizeof(std::complex<double >));
    }

#else
    if (have_local_data) {
        memcpy((void*) data, (void*) multi_array_offset(in, lower, 0, 0),
                local_size_real * sizeof(double));
    }
    fftw_execute(plan);
    if (have_local_data) {
        memcpy((void*) out,
                (void*) (workspace), local_size_complex * sizeof(std::complex<
                        double >));
    }
#endif //USE_FFTW2
}

void
Distributed_fft3d::inv_transform(MArray3dc_ref & in, MArray3d_ref & out)
{
    if (have_local_data) {
        if (in.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible first index offset in input array");
        }
        if (static_cast<int >(in.index_bases()[0] + in.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != get_padded_shape_complex()[1]) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible second dimension of input array");
        }
        if (static_cast<int >(in.shape()[2]) != get_padded_shape_complex()[2]) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible third dimension of input array");
        }
        if (static_cast<int >(out.index_bases()[0]) > lower) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible first index offset in output array");
        }
        if (static_cast<int >(out.index_bases()[0] + out.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible first dimension of output array");
        }
        if (static_cast<int >(out.shape()[1]) != get_padded_shape_real()[1]) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible second dimension of output array");
        }
        if (static_cast<int >(out.shape()[2]) != get_padded_shape_real()[2]) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible third dimension of output array");
        }
    }
#ifdef USE_FFTW2
    if (have_local_data) {
        memcpy( (void*) data, (void*) multi_array_offset(in, lower, 0, 0),
                local_size_complex * sizeof(std::complex<double >));
    }

    rfftwnd_mpi(inv_plan, 1, data, workspace, FFTW_NORMAL_ORDER);

    if (have_local_data) {
        memcpy( (void*) multi_array_offset(out, lower, 0, 0),
                (void*) data, local_size_real * sizeof(double));
    }
#else
    if (have_local_data) {
        memcpy((void*) workspace, (void*) multi_array_offset(in, lower, 0, 0),
                local_size_complex * sizeof(std::complex<double >));
    }
    fftw_execute(inv_plan);
    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0, 0), (void*) data,
                local_size_real * sizeof(double));
    }
#endif //USE_FFTW2
}

void
Distributed_fft3d::inv_transform(fftw_complex * in, MArray3d_ref & out)
{
    if (have_local_data) {
        if (static_cast<int >(out.index_bases()[0]) > lower) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible first index offset in output array");
        }
        if (static_cast<int >(out.index_bases()[0] + out.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible first dimension of output array");
        }
        if (static_cast<int >(out.shape()[1]) != get_padded_shape_real()[1]) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible second dimension of output array");
        }
        if (static_cast<int >(out.shape()[2]) != get_padded_shape_real()[2]) {
            throw std::runtime_error(
                    "Distributed_fft3d::inv_transform found an incompatible third dimension of output array");
        }
    }
#ifdef USE_FFTW2
    if (have_local_data) {
        memcpy( (void*) data, (void*) in,
                local_size_complex * sizeof(std::complex<double >));
    }

    rfftwnd_mpi(inv_plan, 1, data, workspace, FFTW_NORMAL_ORDER);

    if (have_local_data) {
        memcpy( (void*) multi_array_offset(out, lower, 0, 0),
                (void*) data, local_size_real * sizeof(double));
    }
#else
    if (have_local_data) {
        memcpy((void*) workspace, (void*) in,
                local_size_complex * sizeof(std::complex<double >));
    }
    fftw_execute(inv_plan);
    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0, 0), (void*) data,
                local_size_real * sizeof(double));
    }
#endif //USE_FFTW2

}

double
Distributed_fft3d::get_roundtrip_normalization() const
{
    return 1.0 / (shape[0] * shape[1] * shape[2]);
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
