#include "fftw_rectangular_helper.h"
#include "synergia/utils/simple_timer.h"

Fftw_rectangular_helper::Fftw_rectangular_helper(
        std::vector<int > const & grid_shape, Commxx_sptr comm_f_sptr) :
        grid_shape(grid_shape), comm_f_sptr(comm_f_sptr)
{
    if (grid_shape.size() !=3) throw std::runtime_error(
                "Fftw_rectangular_helper: grid_shape does not have 3 elements");
    fftw_mpi_init();
    ptrdiff_t ndim[]={grid_shape[0],grid_shape[1]};


    fftw_local_size = fftw_mpi_local_size_many(2, ndim, grid_shape[2],FFTW_MPI_DEFAULT_BLOCK, comm_f_sptr->get(), &local_nx, &local_x_start);
    if (local_nx == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
    }


    data=(double *) fftw_malloc(sizeof(double) * fftw_local_size);

    fftw_r2r_kind kind_direct[]={FFTW_RODFT10, FFTW_RODFT10};
    plan_r2r_direct=fftw_mpi_plan_many_r2r(2, ndim, grid_shape[2], FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                       data, data, comm_f_sptr->get(), kind_direct, FFTW_EXHAUSTIVE);


   fftw_r2r_kind kind_inv[]={FFTW_RODFT01, FFTW_RODFT01};
   plan_r2r_inv=fftw_mpi_plan_many_r2r(2, ndim, grid_shape[2], FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                       data, data, comm_f_sptr->get(), kind_inv, FFTW_EXHAUSTIVE);

}



 void
Fftw_rectangular_helper::reset_comm_f(Commxx_sptr comm_sptr)
{
    fftw_destroy_plan(plan_r2r_direct);
    fftw_destroy_plan(plan_r2r_inv);
    fftw_free(data);

    ptrdiff_t ndim[]={grid_shape[0],grid_shape[1]};
    fftw_local_size = fftw_mpi_local_size_many(2, ndim, grid_shape[2],FFTW_MPI_DEFAULT_BLOCK, comm_sptr->get(), &local_nx, &local_x_start);
    if (local_nx == 0) {
        have_local_data = false;
    } else {
        have_local_data = true;
    }

    data=(double *) fftw_malloc(sizeof(double) * fftw_local_size);

    fftw_r2r_kind kind_direct[]={FFTW_RODFT10, FFTW_RODFT10};
    plan_r2r_direct=fftw_mpi_plan_many_r2r(2, ndim, grid_shape[2], FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                       data, data, comm_sptr->get(), kind_direct, FFTW_EXHAUSTIVE);


    fftw_r2r_kind kind_inv[]={FFTW_RODFT01, FFTW_RODFT01};
    plan_r2r_inv=fftw_mpi_plan_many_r2r(2, ndim, grid_shape[2], FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,
                                       data, data, comm_sptr->get(), kind_inv, FFTW_EXHAUSTIVE);
}




Fftw_rectangular_helper::~Fftw_rectangular_helper()
{
    fftw_destroy_plan(plan_r2r_direct);
    fftw_destroy_plan(plan_r2r_inv);
    fftw_free(data);
}

ptrdiff_t
Fftw_rectangular_helper::get_local_nx() const
{
    return local_nx;
}

ptrdiff_t
Fftw_rectangular_helper::get_local_x_start() const
{
    return local_x_start;
}

void
Fftw_rectangular_helper::transform(MArray3d_ref & in, MArray3d_ref & out){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (have_local_data) {
        if (static_cast<int >(in.shape()[0]) != local_nx) {
            throw std::runtime_error(
                    "Fftw_rectangular_helper::transformfound an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != grid_shape[1]) {
            throw std::runtime_error(
                    "Fftw_rectangular_helper::transformfound an incompatible second dimension of input array");
        }
        if (static_cast<int >(in.shape()[2]) != grid_shape[2]) {
            throw std::runtime_error(
                    "Fftw_rectangular_helper::transform found an incompatible third dimension of input array");
        }
    }


 //   std::cout<<"local_nx*grid_shape[1]*grid_shape[2]="<<local_nx*grid_shape[1]*grid_shape[2]<<" fftw size="<<fftw_local_size<<std::endl;
    if (have_local_data)  memcpy((void*) data, (void*) in.origin(), local_nx*grid_shape[1]*grid_shape[2]* sizeof(double));
  //  double t;
  //  t = simple_timer_current();
    fftw_execute(plan_r2r_direct);
 //   t = simple_timer_show(t, "fftw_rectangular_helper: transform execute");
    if (have_local_data)  memcpy((void*) out.origin(),( void*) data, local_nx*grid_shape[1]*grid_shape[2]* sizeof(double));
}

void
Fftw_rectangular_helper::inv_transform(MArray3d_ref & in, MArray3d_ref & out){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (have_local_data) {
        if (static_cast<int >(in.shape()[0]) != local_nx) {
            throw std::runtime_error(
                    "Fftw_rectangular_helper::transformfound an incompatible first dimension of input array");
        }
        if (static_cast<int >(in.shape()[1]) != grid_shape[1]) {
            throw std::runtime_error(
                    "Fftw_rectangular_helper::transformfound an incompatible second dimension of input array");
        }
        if (static_cast<int >(in.shape()[2]) != grid_shape[2]) {
            throw std::runtime_error(
                    "Fftw_rectangular_helper::transform found an incompatible third dimension of input array");
        }
    }

    if (have_local_data)
    memcpy((void*) data, (void*) in.origin(), local_nx*grid_shape[1]*grid_shape[2]* sizeof(double));
    fftw_execute(plan_r2r_inv);
    if (have_local_data)
    memcpy((void*) out.origin(),( void*) data, local_nx*grid_shape[1]*grid_shape[2]* sizeof(double));

}
