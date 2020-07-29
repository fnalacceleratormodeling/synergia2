#include <cstring>
#include <stdexcept>
#include "distributed_fft3d_rect_cuda.h"

Distributed_fft3d_rect::Distributed_fft3d_rect()
    : shape()
    , comm(MPI_COMM_NULL)
    , plan()
    , invplan()

    , lower(0)
    , nz(0)
{
}

void Distributed_fft3d_rect::construct(std::array<int, 3> const & new_shape, MPI_Comm new_comm)
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
                "Distributed_fft3d_rect: number of processor must be 1 "
                "for CUDA implementation" );
    }

    shape = new_shape;
    comm = new_comm;

    lower = 0;
    nz = shape[2];

    auto res = cufftPlan3d(&plan, shape[2], shape[1], shape[0], CUFFT_D2Z);
    res = cufftPlan3d(&invplan, shape[2], shape[1], shape[0], CUFFT_Z2D);
}

void
Distributed_fft3d_rect::transform(karray1d_dev& in, karray1d_dev& out)
{
    auto res = cufftExecD2Z( plan, 
            (cufftDoubleReal*)in.data(),
            (cufftDoubleComplex*)out.data() );
}

void
Distributed_fft3d_rect::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
    cufftExecZ2D( invplan, 
            (cufftDoubleComplex*)in.data(),
            (cufftDoubleReal*)out.data() );
}

double
Distributed_fft3d_rect::get_roundtrip_normalization() const
{
    return 1.0 / (shape[0] * shape[1] * shape[2]);
}

Distributed_fft3d_rect::~Distributed_fft3d_rect()
{
    cufftDestroy(plan);
    cufftDestroy(invplan);
}

