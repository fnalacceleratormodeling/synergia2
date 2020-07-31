#include <cstring>
#include <stdexcept>
#include "distributed_fft3d_rect_cuda.h"
#include "synergia/foundation/math_constants.h"

namespace 
{
    void print(karray1d_dev const& arr)
    {
        karray1d harr = Kokkos::create_mirror_view(arr);
        Kokkos::deep_copy(harr, arr);

        std::cout << harr.label() << "\n  ";

        for(int i=0; i<harr.extent(0); ++i)
            std::cout << harr(i) << ", ";

        std::cout << "\n";
    }

    struct alg_zeroer
    {
        karray1d_dev arr;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { arr(i) = 0.0; }
    };

    void zero(karray1d_dev const& arr)
    {
        alg_zeroer alg{arr};
        Kokkos::parallel_for(arr.extent(0), alg);
    }

    struct alg_pad_x
    {
        // in: (x, y, z) of (gx, gy, gz) double
        // out: (z, y, x) of (gz, gy, gx*2) complex
        karray1d_dev in;
        karray1d_dev out;

        int dgx, gy, gz;
        double igygz;
        double igz;

        alg_pad_x(
                karray1d_dev const& in,
                karray1d_dev const& out,
                std::array<int, 3> const& g )
            : in(in), out(out)
            , dgx(2*g[0]), gy(g[1]), gz(g[2])
            , igygz(1.0/(gy*gz))
            , igz(1.0/gz)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int ix = i * igygz;
            int iy = (i - ix*gy*gz) * igz;
            int iz = i - ix*gy*gz - iy*gz;

            // { 1, 2, 3 } -> 
            // { 1, 0, 2, 0, 3, 0, -3, 0, -2, 0, -1, 0 }
            out((iz*gy*dgx + iy*dgx + ix)*2) = in(i);
            out((iz*gy*dgx + iy*dgx + dgx-1-ix)*2) = -in(i);
        }
    };

    struct alg_pad_y
    {
        // in: (z, y, x) of (gz, gy, gx*2) complex
        // out: (z, x, y) of (gz, gx, gy*2) complex
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy;
        int dgx, dgy;
        double igxgy;
        double igx;

        alg_pad_y(
                karray1d_dev const& in,
                karray1d_dev const& out,
                std::array<int, 3> const& g )
            : in(in), out(out)
            , gx(g[0]), gy(g[1])
            , dgx(2*g[0]), dgy(2*g[1])
            , igxgy(1.0/(gy*gx))
            , igx(1.0/gx)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz*gy*gx) * igx;
            int ix = i - iz*gy*gx - iy*gx;

            // { 1, 2, 3 } -> 
            // { 1, 0, 2, 0, 3, 0, -3, 0, -2, 0, -1, 0 }
            out((iz*gx*dgy + ix*dgy + iy)*2) 
                = -in((iz*gy*dgx + iy*dgx + ix+1)*2+1);

            out((iz*gx*dgy + ix*dgy + dgy-1-iy)*2) 
                = in((iz*gy*dgx + iy*dgx + ix+1)*2+1);
        }
    };

    struct alg_pad_z
    {
        // in: (z, x, y) of (gz, gx, gy*2) complex
        // out: (x, y, z) of (gx, gy, gz) double
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy, gz;
        int dgy;
        double igxgy;
        double igx;

        alg_pad_z(
                karray1d_dev const& in,
                karray1d_dev const& out,
                std::array<int, 3> const& g )
            : in(in), out(out)
            , gx(g[0]), gy(g[1]), gz(g[2])
            , dgy(g[1]*2)
            , igxgy(1.0/(gx*gy))
            , igx(1.0/gx)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz*gy*gx) * igx;
            int ix = i - iz*gy*gx - iy*gx;

            out(ix*gy*gz + iy*gz + iz)
                = -in((iz*gx*dgy + ix*dgy + iy+1)*2 + 1);
        }
    };


    struct alg_shift
    {
        karray1d_dev in;
        int gx;
        double igx;

        alg_shift(karray1d_dev const& in, int gx)
            : in(in), gx(gx), igx(1.0/gx)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            const int idx = i * igx;
            const int ix = i - idx*gx;

            double dk = (ix+1) * mconstants::pi * igx * 0.5;

            double imag = 
                - in(idx*gx*4 + (ix+1)*2 + 0) * sin(dk)
                + in(idx*gx*4 + (ix+1)*2 + 1) * cos(dk);

            in(idx*gx*4 + (ix+1)*2 + 1) = imag;
        }
    };

}

Distributed_fft3d_rect::Distributed_fft3d_rect()
    : shape()
    , comm(MPI_COMM_NULL)
    , datax()
    , datay()
    , plan_x()
    , plan_y()
    , plan_z()
    , lower(0)
    , nx(0)
{
}

void Distributed_fft3d_rect::construct(std::array<int, 3> const & new_shape, MPI_Comm new_comm)
{
    cufftDestroy(plan_x);
    cufftDestroy(plan_y);
    cufftDestroy(plan_z);

    //cufftDestroy(invplan);

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
    nx = shape[2];

    // (x*2, y, z) of a complex grid for batched DST along x
    // (x, y*2, z) of a complex grid for batched DST along y
    datax = karray1d_dev("datax", shape[0]*4*shape[1]*shape[2]);
    datay = karray1d_dev("datay", shape[0]*shape[1]*4*shape[2]);

    // plan x
    int rank = 1;
    int istride = 1;
    int ostride = 1;

    int n_x[] = {shape[0]*2};
    int inembed_x[] = {shape[0]*2};
    int onembed_x[] = {shape[0]*2};
    int idist_x = shape[0]*2;
    int odist_x = shape[0]*2;

    cufftPlanMany(&plan_x, rank, n_x,
            inembed_x, istride, idist_x,
            onembed_x, ostride, odist_x,
            CUFFT_Z2Z, shape[1]*shape[2]);

    // plan y
    int n_y[] = {shape[1]*2};
    int inembed_y[] = {shape[1]*2};
    int onembed_y[] = {shape[1]*2};
    int idist_y = shape[1]*2;
    int odist_y = shape[1]*2;

    cufftPlanMany(&plan_y, rank, n_y,
            inembed_y, istride, idist_y,
            onembed_y, ostride, odist_y,
            CUFFT_Z2Z, shape[0]*shape[2]);

    // plan z
    int n_z[] = {shape[2]};
    int inembed_z[] = {shape[2]};
    int onembed_z[] = {shape[2]};
    int idist_z = shape[2];
    int odist_z = shape[2]/2+1;

    cufftPlanMany(&plan_z, rank, n_z,
            inembed_z, istride, idist_z,
            onembed_z, ostride, odist_z,
            CUFFT_D2Z, shape[0]*shape[1]);
}

void
Distributed_fft3d_rect::transform(karray1d_dev& in, karray1d_dev& out)
{
    const int gx = shape[0];
    const int gy = shape[1];
    const int gz = shape[2];

    zero(datax);
    zero(datay);

    alg_pad_x padx(in, datax, shape);
    Kokkos::parallel_for(gx*gy*gz, padx);

    auto res = cufftExecZ2Z(plan_x,
            (cufftDoubleComplex*)datax.data(),
            (cufftDoubleComplex*)datax.data(),
            CUFFT_FORWARD);

    alg_shift shiftx(datax, gx);
    Kokkos::parallel_for(gx*gy*gz, shiftx);

    alg_pad_y pady(datax, datay, shape);
    Kokkos::parallel_for(gx*gy*gz, pady);

    res = cufftExecZ2Z(plan_y,
            (cufftDoubleComplex*)datay.data(),
            (cufftDoubleComplex*)datay.data(),
            CUFFT_FORWARD);

    alg_shift shifty(datay, gy);
    Kokkos::parallel_for(gx*gy*gz, shifty);

    alg_pad_z padz(datay, datax, shape);
    Kokkos::parallel_for(gx*gy*gz, padz);

    res = cufftExecD2Z(plan_z,
            (cufftDoubleReal*)datax.data(),
            (cufftDoubleComplex*)out.data() );
}

void
Distributed_fft3d_rect::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
#if 0
    cufftExecZ2D( invplan, 
            (cufftDoubleComplex*)in.data(),
            (cufftDoubleReal*)out.data() );
#endif
}

double
Distributed_fft3d_rect::get_roundtrip_normalization() const
{
    return 1.0 / (shape[0] * shape[1] * shape[2]);
}

Distributed_fft3d_rect::~Distributed_fft3d_rect()
{
    cufftDestroy(plan_x);
    cufftDestroy(plan_y);
    cufftDestroy(plan_z);

    //cufftDestroy(invplan);
}

