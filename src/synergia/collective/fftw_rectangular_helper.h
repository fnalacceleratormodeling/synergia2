#ifndef FFTW_RECTANGULAR_HELPER_H_
#define FFTW_RECTANGULAR_HELPER_H_
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <boost/shared_ptr.hpp>
#include "synergia/utils/commxx.h"
#include <vector>
#include <stdexcept>
#include <string.h> // needed for memcpy on hopper <string> won't work
#include "synergia/utils/multi_array_typedefs.h"

class Fftw_rectangular_helper
{   std::vector<int > grid_shape;
    fftw_plan  plan_r2r_direct, plan_r2r_inv;
    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size;
    double *data;
    Commxx_sptr comm_f_sptr;
    bool have_local_data;


public:

    Fftw_rectangular_helper( std::vector<int >  const & grid_shape, Commxx_sptr comm_f_sptr);

    void
    reset_comm_f(Commxx_sptr comm_sptr);

    ptrdiff_t
    get_local_nx() const;

    ptrdiff_t
    get_local_x_start() const;

    void
    transform(MArray3d_ref & in, MArray3d_ref & out);

    void
    inv_transform(MArray3d_ref & in, MArray3d_ref & out);

    ~Fftw_rectangular_helper();
};

typedef boost::shared_ptr<Fftw_rectangular_helper> Fftw_rectangular_helper_sptr; // syndoc:include
#endif /* FFTW_RECTANGULAR_HELPER_H_ */
