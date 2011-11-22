#include "distributed_rectangular_grid.h"
#include "synergia/utils/multi_array_offsets.h"

void
Distributed_rectangular_grid::construct_hockney(int lower, int upper,
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
    if (lower == upper) {
        lower_guard = lower;
        upper_guard = upper;
    }
    grid_points_sptr
            = boost::shared_ptr<MArray3d >(
                    new MArray3d(boost::extents[extent_range(lower_guard,
                            upper_guard)][array_shape[1]][array_shape[2]]));
    grid_points_2dc_sptr
            = boost::shared_ptr<MArray2dc >(
                    new MArray2dc(boost::extents[extent_range(lower_guard,
                            upper_guard)][array_shape[1]]));
    grid_points_1d_sptr
            = boost::shared_ptr<MArray1d >(
                    new MArray1d(boost::extents[array_shape[2]]));
    normalization = 1.0;
}

void
Distributed_rectangular_grid::construct_rectangular(int lower, int upper,
        std::vector<int > const & array_shape)
{
    std::vector<int > grid_shape(domain_sptr->get_grid_shape());
    this->lower = lower;
    this->upper = upper;
    if (lower == 0)  {
        lower_guard = 0;
    } else {
        lower_guard = lower - 1;
    }
    if (upper == grid_shape[0]) {
        upper_guard = grid_shape[0];
    } else {
        upper_guard = upper + 1;
    }
    if (lower == upper) {
        lower_guard = lower;
        upper_guard = upper;
    }


    grid_points_sptr
            = boost::shared_ptr<MArray3d >(
                    new MArray3d(boost::extents[array_shape[0]][array_shape[1]][array_shape[2]]));  
 
     grid_points_2dc_sptr
              = boost::shared_ptr<MArray2dc >(
                      new MArray2dc(boost::extents[array_shape[0]][array_shape[1]]));

    grid_points_1d_sptr
            = boost::shared_ptr<MArray1d >(
                    new MArray1d(boost::extents[array_shape[2]]));
    normalization = 1.0;
}



Distributed_rectangular_grid::Distributed_rectangular_grid(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic, int lower,
        int upper, Commxx const& comm, std::string const solver) :
    comm(comm), uppers(0), lengths(0)
{
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape, periodic));
    if (solver=="hockney") {        
        construct_hockney(lower, upper, domain_sptr->get_grid_shape());
    }
    else if  (solver=="rectangular"){    
         construct_rectangular(lower, upper, domain_sptr->get_grid_shape());
    }    

}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        Rectangular_grid_domain_sptr rectangular_grid_domain_sptr, int lower,
        int upper, Commxx const& comm, std::string const solver) :
    comm(comm), uppers(0), lengths(0)
{
    domain_sptr = rectangular_grid_domain_sptr;
    if (solver=="hockney") {
        construct_hockney(lower, upper, domain_sptr->get_grid_shape());
    }
    else if  (solver=="rectangular"){    
         construct_rectangular(lower, upper, domain_sptr->get_grid_shape());
    }   
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        Rectangular_grid_domain_sptr rectangular_grid_domain_sptr, int lower,
        int upper, std::vector<int > const & padded_shape, Commxx const& comm) :
    comm(comm), uppers(0), lengths(0)
{
    domain_sptr = rectangular_grid_domain_sptr;
    construct_hockney(lower, upper, padded_shape);
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

void
Distributed_rectangular_grid::calculate_uppers_lengths()
{
    if (uppers.empty()) {
        int size = comm.get_size();
        uppers.resize(size);
        lengths.resize(size);
        if (size == 1) {
            uppers[0] = upper;
        } else {
            MPI_Allgather((void*) (&upper), 1, MPI_INT, (void*) (&uppers[0]),
                    1, MPI_INT, comm.get());
        }
        std::vector<int> shape(domain_sptr->get_grid_shape());
        lengths[0] = uppers[0] * shape[1] * shape[2];
        for (int i = 1; i < size; ++i) {
            if (uppers[i] == 0) {
                uppers[i] = uppers[i - 1];
            }
            lengths[i] = (uppers[i] - uppers[i - 1]) * shape[1] * shape[2];
        }
    }
}

std::vector<int > const&
Distributed_rectangular_grid::get_uppers()
{
    calculate_uppers_lengths();
    return uppers;
}

std::vector<int > const&
Distributed_rectangular_grid::get_lengths()
{
    calculate_uppers_lengths();
    return lengths;
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

MArray2dc_ref const&
Distributed_rectangular_grid::get_grid_points_2dc() const
{
    return *grid_points_2dc_sptr;
}

MArray2dc_ref &
Distributed_rectangular_grid::get_grid_points_2dc()
{
    return *grid_points_2dc_sptr;
}

MArray1d_ref const&
Distributed_rectangular_grid::get_grid_points_1d() const
{
    return *grid_points_1d_sptr;
}

MArray1d_ref &
Distributed_rectangular_grid::get_grid_points_1d()
{
    return *grid_points_1d_sptr;
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

Commxx const &
Distributed_rectangular_grid::get_comm() const
{
return comm;
}

void
Distributed_rectangular_grid::fill_guards()
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
    recv_buffer = (void*) multi_array_offset(*grid_points_sptr, lower_guard, 0,
            0);
    send_buffer
            = (void*) multi_array_offset(*grid_points_sptr, upper - 1, 0, 0);
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
    if (!domain_sptr->is_periodic()) {
        if (upper >= domain_sptr->get_grid_shape()[0]) {
            send = false;
        }
        if ((lower >= domain_sptr->get_grid_shape()[0])) {
            receive = false;
        }
    }
    if (send) {
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
                comm.get());
    }
    if (receive) {
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, sender, sender,
                comm.get(), &status);
    }

    //send to the left
    recv_buffer = (void*) multi_array_offset(*grid_points_sptr,
            upper_guard - 1, 0, 0);
    send_buffer = (void*) multi_array_offset(*grid_points_sptr, lower, 0, 0);
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
    if (!domain_sptr->is_periodic()) {
        if (lower >= domain_sptr->get_grid_shape()[0]) {
            send = false;
        }
        if ((upper >= domain_sptr->get_grid_shape()[0]) && !(rank == size - 1)) {
            receive = false;
        }
    }
    if (send) {
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
                comm.get());
    }
    if (receive) {
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, sender, sender,
                comm.get(), &status);
    }
}
