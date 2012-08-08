#include <iostream>
#include <cstring>
#include <stdexcept>
#include "multi_array_offsets.h"
#include "distributed_fft2d.h"
#include "synergia/utils/simple_timer.h"

Distributed_fft2d::Distributed_fft2d(std::vector<int > const & shape,
        Commxx_sptr comm_sptr, int planner_flags,
        std::string const& wisdom_filename) :
        uppers(0), lengths(0), lengths_1d(0), shape(shape), comm_sptr(comm_sptr)
{
    if (comm_sptr->get_size() / 2 >= shape[0] / 2) {
        throw std::runtime_error(
                "Distributed_fft2d: (number of processors)/2 must be <= shape[0]/2");
    }
    //n.b. : we aren't using the wisdom_filename yet
#ifdef USE_FFTW2
    plan = fftwnd_mpi_create_plan(comm_sptr->get(), 2, &shape[0],
            FFTW_FORWARD, planner_flags);
    inv_plan = fftwnd_mpi_create_plan(comm_sptr->get(), 2, &shape[0],
            FFTW_BACKWARD, planner_flags);
    // fftw2 mpi manual says to allocate data and workspace of
    // size total_local_size, not local_nx size.
    // http://www.fftw.org/fftw2_doc/fftw_4.html#SEC47
    int total_local_size;
    int local_nx, local_x_start;
    int local_ny_after_transpose, local_y_start_after_transpose;
    fftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
            &local_ny_after_transpose,
            &local_y_start_after_transpose,
            &total_local_size);

    data = (fftw_complex *) fftw_malloc(total_local_size * sizeof(fftw_complex));
    workspace = (fftw_complex *) fftw_malloc(total_local_size * sizeof(fftw_complex));
    local_size_real = local_nx * shape[1];
    if (local_size_real == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
    }
#else
    fftw_mpi_init();
    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_2d(shape[0], shape[1],
            comm_sptr->get(), &local_nx, &local_x_start);
    local_size_real = local_nx * shape[1];
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
    //data = (fftw_complex *) fftw_alloc_complex(sizeof(fftw_complex) * fftw_local_size);
    //workspace = (fftw_complex *) fftw_alloc_complex(sizeof(fftw_complex) * fftw_local_size);
    //data = fftw_alloc_complex(fftw_local_size);
    //workspace = fftw_alloc_complex(fftw_local_size);
    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftw_local_size);
    workspace = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftw_local_size);

    plan = fftw_mpi_plan_dft_2d(shape[0], shape[1], data, workspace,
            comm_sptr->get(), FFTW_FORWARD, planner_flags);
    inv_plan = fftw_mpi_plan_dft_2d(shape[0], shape[1], workspace, data,
            comm_sptr->get(), FFTW_BACKWARD, planner_flags);
#endif //USE_FFTW2
    lower = local_x_start;
    upper = lower + local_nx;
}

int
Distributed_fft2d::get_data_size() const
{
    return sizeof(data);
}

int
Distributed_fft2d::get_workspace_size() const
{
    return sizeof(workspace);
}

Commxx &
Distributed_fft2d::get_comm()
{
    return *comm_sptr;
}

Commxx_sptr
Distributed_fft2d::get_comm_sptr()
{
    return comm_sptr;
}

int
Distributed_fft2d::get_lower() const
{
    return lower;
}

int
Distributed_fft2d::get_upper() const
{
    return upper;
}

void
Distributed_fft2d::calculate_uppers_lengths()
{
    if (uppers.empty()) {
        int size = comm_sptr->get_size();
        uppers.resize(size);
        lengths.resize(size);
        lengths_1d.resize(size);
        if (size == 1) {
            uppers[0] = upper;
        } else {
            MPI_Allgather((void*) (&upper), 1, MPI_INT, (void*) (&uppers[0]),
                    1, MPI_INT, comm_sptr->get());
        }
        lengths[0] = uppers[0] * shape[1];
        lengths_1d[0] = shape[2];
        for (int i = 1; i < size; ++i) {
            if (uppers[i] == 0) {
                uppers[i] = uppers[i - 1];
            }
            lengths[i] = (uppers[i] - uppers[i - 1]) * shape[1];
            lengths_1d[i] = shape[2];
        }
    }
}

std::vector<int > const&
Distributed_fft2d::get_uppers()
{
    calculate_uppers_lengths();
    return uppers;
}

std::vector<int > const&
Distributed_fft2d::get_lengths()
{
    calculate_uppers_lengths();
    return lengths;
}

std::vector<int > const&
Distributed_fft2d::get_lengths_1d()
{
    calculate_uppers_lengths();
    return lengths_1d;
}

std::vector<int > const&
Distributed_fft2d::get_shape() const
{
    return shape;
}

void
Distributed_fft2d::transform(MArray2dc_ref & in, MArray2dc_ref & out)
{
    //double t;
    //t = simple_timer_current();
    if (have_local_data) {
        if (in.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft2d::transform found an incompatible first index offset in input array");
        }
        if (static_cast<int >(in.index_bases()[0] + in.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft2d::transform found an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != get_shape()[1]) {
            throw std::runtime_error(
                    "Distributed_fft2d::transform found an incompatible second dimension of input array");
        }
        if (out.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft2d::transform found an incompatible first index offset in output array");
        }
        if (static_cast<int >(out.index_bases()[0] + out.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft2d::transform found an incompatible first dimension of output array");
        }
        if (static_cast<int >(out.shape()[1]) != get_shape()[1]) {
            throw std::runtime_error(
                    "Distributed_fft2d::transform found an incompatible second dimension of output array");
        }
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(error-check)");
#ifdef USE_FFTW2
    if (have_local_data) {
        memcpy((void*) data, (void*) multi_array_offset(in, lower, 0),
	       local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(memcpy-in)");
    fftwnd_mpi(plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(fftw_execute)");
    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0),
	       (void*) (data), local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(memcpy-out)");
#else
    if (have_local_data) {
        memcpy((void*) data, (void*) multi_array_offset(in, lower, 0),
                local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(memcpy-in)");
    fftw_execute(plan);
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(fftw_execute)");
    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0),
                (void*) (workspace), local_size_real * sizeof(std::complex<
                        double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-transform(memcpy-out)");
#endif //USE_FFTW2
}

void
Distributed_fft2d::inv_transform(MArray2dc_ref & in, MArray2dc_ref & out)
{
    double t;
    //t = simple_timer_current();
    if (have_local_data) {
        if (in.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft2d::inv_transform found an incompatible first index offset in input array");
        }
        if (static_cast<int >(in.index_bases()[0] + in.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft2d::inv_transform found an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != get_shape()[1]) {
            throw std::runtime_error(
                    "Distributed_fft2d::inv_transform found an incompatible second dimension of input array");
        }
        if (out.index_bases()[0] > lower) {
            throw std::runtime_error(
                    "Distributed_fft2d::inv_transform found an incompatible first index offset in output array");
        }
        if (static_cast<int >(out.index_bases()[0] + out.shape()[0]) < upper) {
            throw std::runtime_error(
                    "Distributed_fft2d::inv_transform found an incompatible first dimension of output array");
        }
        if (static_cast<int >(out.shape()[1]) != get_shape()[1]) {
            throw std::runtime_error(
                    "Distributed_fft2d::inv_transform found an incompatible second dimension of output array");
        }
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(error-check)");
#ifdef USE_FFTW2
    if (have_local_data) {
        memcpy((void*) data, (void*) multi_array_offset(in, lower, 0),
	       local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(memcpy-in)");
    fftwnd_mpi(inv_plan, 1, data, workspace, FFTW_NORMAL_ORDER);
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(fftw_excute)");

    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0), (void*) data,
	       local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(memcpy-out)");
#else
    if (have_local_data) {
        memcpy((void*) workspace, (void*) multi_array_offset(in, lower, 0),
                local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(memcpy-in)");
    fftw_execute(inv_plan);
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(fftw_excute)");

    if (have_local_data) {
        memcpy((void*) multi_array_offset(out, lower, 0), (void*) data,
                local_size_real * sizeof(std::complex<double >));
    }
    //t = simple_timer_show(t, "sc-distributed_fft2d-inv_transform(memcpy-out)");
#endif //USE_FFTW2
}

double
Distributed_fft2d::get_roundtrip_normalization() const
{
    return 1.0 / (shape[0] * shape[1] );
}

Distributed_fft2d::~Distributed_fft2d()
{
#ifdef USE_FFTW2
    fftwnd_mpi_destroy_plan(plan);
    fftwnd_mpi_destroy_plan(inv_plan);
    fftw_free(data);
    fftw_free(workspace);
#else
    fftw_destroy_plan(plan);
    fftw_destroy_plan(inv_plan);
    fftw_free(data);
    fftw_free(workspace);
#endif //USE_FFTW2
}
