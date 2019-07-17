#include <cstring>
#include <stdexcept>
#include "distributed_fft2d.h"

Distributed_fft2d::Distributed_fft2d(std::array<int, 3> const & shape)
    : shape(shape)
    , fft(shape)

    , plan()
    , inv_plan()
    , data(nullptr)
    , workspace(nullptr)
    , lower(0)
    , upper(0)
    , uppers(0)
    , lengths(0)
    , lengths_1d(0)
    , local_size_real(0)
    , have_local_data(false)
{
#if 0
    if (comm_sptr->get_size() / 2 >= shape[0] / 2) {
        throw std::runtime_error(
                "Distributed_fft2d: (number of processors)/2 must be <= shape[0]/2");
    }
#endif

    //n.b. : we aren't using the wisdom_filename yet
    fftw_mpi_init();

    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_2d(
            shape[0], shape[1], MPI_COMM_WORLD,
            &local_nx, &local_x_start);

    local_size_real = local_nx * shape[1];
    have_local_data = (local_nx > 0);

    // MEDIUM HACK. fftw often (always?) returns a local size that is
    // impossibly small. Adjust to the smallest possible size.
    if (fftw_local_size < local_size_real)
        fftw_local_size = local_size_real;

    data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftw_local_size);
    workspace = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftw_local_size);

    plan = fftw_mpi_plan_dft_2d( 
            shape[0], shape[1], data, workspace,
            MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE );

    inv_plan = fftw_mpi_plan_dft_2d(
            shape[0], shape[1], workspace, data,
            MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE );

    lower = local_x_start;
    upper = lower + local_nx;
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
#if 0
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
#endif
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

std::array<int, 3> const&
Distributed_fft2d::get_shape() const
{
    return shape;
}

void
Distributed_fft2d::transform(karray1d_dev & in, karray1d_dev & out)
{
    fft.transform(in, out);

#if 0
    if (have_local_data) 
    {
        memcpy( (void*)data, 
                (void*)&in(lower*shape[1]*2),
                local_size_real * sizeof(double) * 2);
    }

    fftw_execute(plan);

    if (have_local_data) 
    {
        memcpy( (void*)&out(lower*shape[1]*2),
                (void*)(workspace), 
                local_size_real * sizeof(double) * 2);
    }
#endif
}

void
Distributed_fft2d::inv_transform(karray1d_dev & in, karray1d_dev & out)
{
    fft.inv_transform(in, out);

#if 0
    if (have_local_data) 
    {
        memcpy( (void*)workspace, 
                (void*)&in(lower*shape[1]*2),
                local_size_real * sizeof(double) * 2);
    }

    fftw_execute(inv_plan);

    if (have_local_data) 
    {
        memcpy( (void*)&out(lower*shape[1]*2),
                (void*)data,
                local_size_real * sizeof(double) * 2);
    }
#endif
}

double
Distributed_fft2d::get_roundtrip_normalization() const
{
    return 1.0 / (shape[0] * shape[1] );
}

Distributed_fft2d::~Distributed_fft2d()
{
#if 0
    fftw_destroy_plan(plan);
    fftw_destroy_plan(inv_plan);
    fftw_free(data);
    fftw_free(workspace);
#endif
}

