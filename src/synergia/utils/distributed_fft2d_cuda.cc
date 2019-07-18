#include <cstring>
#include <stdexcept>
#include "distributed_fft2d.h"

Distributed_fft2d::Distributed_fft2d(std::array<int, 3> const & shape, MPI_Comm comm)
    : shape(shape)
    , comm(comm)
    , plan()

    , lower(0)
    , upper(shape[0])

    , uppers(1)
    , lengths(1)
    , lengths_1d(1)
{
    int comm_size = 0;
    MPI_Comm_size(comm, &comm_size);

    if (comm_size != 1)
    {
        throw std::runtime_error(
                "Distributed_fft2d: number of processor must be 1 "
                "for CUDA implementation" );
    }

    uppers[0] = shape[0];
    lengths[0] = shape[0] * shape[1];
    lengths_1d[0] = shape[2];

    cufftPlan2d(&plan, shape[0], shape[1], CUFFT_Z2Z);
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

std::vector<int > const&
Distributed_fft2d::get_uppers()
{
    return uppers;
}

std::vector<int > const&
Distributed_fft2d::get_lengths()
{
    return lengths;
}

std::vector<int > const&
Distributed_fft2d::get_lengths_1d()
{
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
    cufftExecZ2Z( plan, 
            (cufftDoubleComplex*)in.data(),
            (cufftDoubleComplex*)out.data(),
            CUFFT_FORWARD );
}

void
Distributed_fft2d::inv_transform(karray1d_dev & in, karray1d_dev & out)
{
    cufftExecZ2Z( plan, 
            (cufftDoubleComplex*)in.data(),
            (cufftDoubleComplex*)out.data(),
            CUFFT_INVERSE );
}

double
Distributed_fft2d::get_roundtrip_normalization() const
{
    return 1.0 / (shape[0] * shape[1] );
}

Distributed_fft2d::~Distributed_fft2d()
{
    cufftDestroy(plan);
}

