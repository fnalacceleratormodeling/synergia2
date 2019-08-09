#include <cstring>
#include <stdexcept>
#include "distributed_fft2d_cuda.h"

Distributed_fft2d::Distributed_fft2d()
    : shape()
    , comm(MPI_COMM_NULL)
    , plan()

    , lower(0)
    , nx(0)
{
}

void Distributed_fft2d::construct(std::array<int, 3> const & new_shape, MPI_Comm new_comm)
{
    cufftDestroy(plan);

    if (new_comm == MPI_COMM_NULL)
    {
        comm = MPI_COMM_NULL;
        return;
    }

    int comm_size = 0;
    MPI_Comm_size(new_comm, &comm_size);

    if (comm_size != 1)
    {
        throw std::runtime_error(
                "Distributed_fft2d: number of processor must be 1 "
                "for CUDA implementation" );
    }

    shape = new_shape;
    comm = new_comm;

    lower = 0;
    nx = shape[0];

    cufftPlan2d(&plan, shape[0], shape[1], CUFFT_Z2Z);
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

