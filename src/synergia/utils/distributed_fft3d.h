#ifndef DISTRIBUTED_FFT3D_H_
#define DISTRIBUTED_FFT3D_H_

#ifdef USE_FFTW2
#include <rfftw_mpi.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif //USE_FFTW2
#include <vector>
#include <string>
#include "boost/shared_ptr.hpp"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

class Distributed_fft3d
{
private:
#ifdef USE_FFTW2
    rfftwnd_mpi_plan plan, inv_plan;
    fftw_real *data;
    fftw_real *workspace;
#else
    fftw_plan plan, inv_plan;
    double *data;
    fftw_complex *workspace;
    fftw_complex *workspace1;
    fftw_complex *workspace2;
#endif //USE_FFTW2
    int lower, upper;
    std::vector<int> uppers, lengths;
    int local_size_real, local_size_complex;
    std::vector<int > shape;
    bool have_local_data;

    int nt;   // number of threads
    int np;   // number of mpi procs
    int rank; // mpi rank

    int n0, n1, n2, n2r, n2c;
    fftw_plan fp1, fp2, fp3;
    fftw_plan bp1, bp2, bp3;

    // complex type (contiguous of two doubles)
    MPI_Datatype mpi_cplx;

    // data types for the all to all BEFORE 3rd dim fft
    MPI_Datatype mpi_fs11, mpi_fs12;
    MPI_Datatype mpi_fr11, mpi_fr12, mpi_fr13, mpi_fr14;

    // data types for the all to all AFTER 3rd dim fft
    MPI_Datatype mpi_fs21, mpi_fs22, mpi_fs23;
    MPI_Datatype mpi_fr21, mpi_fr22;

    Commxx_sptr comm_sptr;
    void
    calculate_uppers_lengths();
public:
    Distributed_fft3d(std::vector<int > const& shape, Commxx_sptr comm_sptr,
            int planner_flags = FFTW_ESTIMATE,
            std::string const& wisdom_filename = "");
    Commxx_sptr
    get_comm_sptr();
    Commxx &
    get_comm();
    int
    get_nt() const { return nt; }
    int
    get_np() const { return np; }
    int
    get_rank() const { return rank; }
    int
    get_lower() const;
    int
    get_upper() const;
    std::vector<int > const&
    get_uppers();
    std::vector<int > const&
    get_lengths();
    std::vector<int > const&
    get_shape() const;
    std::vector<int >
    get_padded_shape_real() const;
    std::vector<int >
    get_padded_shape_complex() const;
    void
    transform(MArray3d_ref & in, MArray3dc_ref & out);
    void
    transform(MArray3d_ref & in, fftw_complex * out);
    void
    inv_transform(MArray3dc_ref & in, MArray3d_ref & out);
    void
    inv_transform(fftw_complex * in, MArray3d_ref & out);
    double
    get_roundtrip_normalization() const;
    ~Distributed_fft3d();
};

typedef boost::shared_ptr<Distributed_fft3d > Distributed_fft3d_sptr; // syndoc:include
#endif /* DISTRIBUTED_FFT3D_H_ */
