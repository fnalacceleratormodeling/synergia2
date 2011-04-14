#include "distributed_rectangular_grid.h"

void
Distributed_rectangular_grid::construct(int lower, int upper,
        std::vector<int > const & array_shape)
{
    std::vector<int > grid_shape(domain_sptr->get_grid_shape());
    this->lower = lower;
    this->upper = upper;
    if ((lower == 0) && (!domain_sptr->is_periodic())) {
        lower_guard = 0;
    } else {
        lower_guard = lower - 1;
    }
    if ((upper == grid_shape[0]) && (!domain_sptr->is_periodic())) {
        upper_guard = grid_shape[0];
    } else {
        upper_guard = upper + 1;
    }
    grid_points_sptr
            = boost::shared_ptr<MArray3d >(
                    new MArray3d(boost::extents[extent_range(lower_guard,
                            upper_guard)][array_shape[1]][array_shape[2]]));
    normalization = 1.0;
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic, int lower,
        int upper, Commxx const& commxx) :
    commxx(commxx)
{
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape, periodic));
    construct(lower, upper, domain_sptr->get_grid_shape());
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        Rectangular_grid_domain_sptr rectangular_grid_domain_sptr, int lower,
        int upper, Commxx const& commxx) :
    commxx(commxx)
{
    domain_sptr = rectangular_grid_domain_sptr;
    construct(lower, upper, domain_sptr->get_grid_shape());
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        Rectangular_grid_domain_sptr rectangular_grid_domain_sptr, int lower,
        int upper, std::vector<int > const & padded_shape, Commxx const& commxx) :
    commxx(commxx)
{
    domain_sptr = rectangular_grid_domain_sptr;
    construct(lower, upper, padded_shape);
}

Rectangular_grid_domain_sptr
Distributed_rectangular_grid::get_domain_sptr() const
{
    return domain_sptr;
}

Rectangular_grid_domain_sptr
Distributed_rectangular_grid::get_domain_sptr()
{
    return domain_sptr;
}

int
Distributed_rectangular_grid::get_lower() const
{
    return lower;
}

int
Distributed_rectangular_grid::get_upper() const
{
    return upper;
}

int
Distributed_rectangular_grid::get_lower_guard() const
{
    return lower_guard;
}

int
Distributed_rectangular_grid::get_upper_guard() const
{
    return upper_guard;
}

MArray3d_ref const&
Distributed_rectangular_grid::get_grid_points() const
{
    return *grid_points_sptr;
}

MArray3d_ref &
Distributed_rectangular_grid::get_grid_points()
{
    return *grid_points_sptr;
}

void
Distributed_rectangular_grid::set_normalization(double val)
{
    normalization = val;
}

double
Distributed_rectangular_grid::get_normalization() const
{
    return normalization;
}

void
Distributed_rectangular_grid::fill_guards(Commxx const & comm)
{
    int rank = comm.get_rank();
    int size = comm.get_size();
    if (size == 1) {
        if (domain_sptr->is_periodic()) {
            int max0 = domain_sptr->get_grid_shape()[0];
            int max1 = domain_sptr->get_grid_shape()[1];
            int max2 = domain_sptr->get_grid_shape()[2];
            for (int j = 0; j < max1; ++j) {
                for (int k = 0; k < max2; ++k) {
                    (*grid_points_sptr)[-1][j][k] = (*grid_points_sptr)[max0
                            - 1][j][k];
                }
            }
            for (int j = 0; j < max1; ++j) {
                for (int k = 0; k < max2; ++k) {
                    (*grid_points_sptr)[max0][j][k]
                            = (*grid_points_sptr)[0][j][k];
                }
            }
            return;
        } else {
            return;
        }
    }

    MPI_Status status;
    void *recv_buffer, *send_buffer;
    // send to the right
    recv_buffer = (void*) (grid_points_sptr->origin() + lower_guard
            * grid_points_sptr->strides()[0]);
    send_buffer = (void*) (grid_points_sptr->origin() + (upper - 1)
            * grid_points_sptr->strides()[0]);
    size_t message_size = grid_points_sptr->shape()[1]
            * grid_points_sptr->shape()[2];
    int sender = rank - 1;
    bool send = true;
    int receiver = rank + 1;
    bool receive = true;
    if (rank == size - 1) {
        if (domain_sptr->is_periodic()) {
            receiver = 0;
        } else {
            send = false;
        }
    }
    if (rank == 0) {
        if (domain_sptr->is_periodic()) {
            sender = size - 1;
        } else {
            receive = false;
        }
    }

    if (send) {
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
                commxx.get());
    }
    if (receive) {
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, sender, sender,
                commxx.get(), &status);
    }

    //send to the left
    recv_buffer = (void*) (grid_points_sptr->origin() + (upper_guard - 1)
            * grid_points_sptr->strides()[0]);
    send_buffer = (void*) (grid_points_sptr->origin() + lower
            * grid_points_sptr->strides()[0]);
    sender = rank + 1;
    send = true;
    receiver = rank - 1;
    receive = true;
    if (rank == size - 1) {
        if (domain_sptr->is_periodic()) {
            sender = 0;
        } else {
            receive = false;
        }
    }
    if (rank == 0) {
        if (domain_sptr->is_periodic()) {
            receiver = size - 1;
        } else {
            send = false;
        }
    }

    if (send) {
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
                commxx.get());
    }
    if (receive) {
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, sender, sender,
                commxx.get(), &status);
    }

}
