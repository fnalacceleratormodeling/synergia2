#ifndef DISTRIBUTED_FFT2D_H_
#define DISTRIBUTED_FFT2D_H_

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

class Distributed_fft2d
{
private:
#ifdef USE_FFTW2
    fftwnd_mpi_plan plan, inv_plan;
    fftw_complex *data;
    fftw_complex *workspace;
#else
    fftw_plan plan, inv_plan;
    fftw_complex *data;
    fftw_complex *workspace;
#endif //USE_FFTW2
    int lower, upper;
    std::vector<int> uppers, lengths, lengths_1d;
    int local_size_real;
    std::vector<int > shape;
    bool have_local_data;
    Commxx comm;
    void
    calculate_uppers_lengths();
public:
    Distributed_fft2d(std::vector<int > const& shape, Commxx const& comm,
            int planner_flags = FFTW_ESTIMATE,
            std::string const& wisdom_filename = "");

    int
    get_data_size() const;
    int
    get_workspace_size() const;


    Commxx &
    get_comm();
    int
    get_lower() const;
    int
    get_upper() const;
    std::vector<int > const&
    get_uppers();
    std::vector<int > const&
    get_lengths();
    std::vector<int > const&
    get_lengths_1d();
    std::vector<int > const&
    get_shape() const;
    void
    transform(MArray2dc_ref & in, MArray2dc_ref & out);
    void
    inv_transform(MArray2dc_ref & in, MArray2dc_ref & out);
    double
    get_roundtrip_normalization() const;
    ~Distributed_fft2d();
};

typedef boost::shared_ptr<Distributed_fft2d > Distributed_fft2d_sptr;
#endif /* DISTRIBUTED_FFT2D_H_ */
