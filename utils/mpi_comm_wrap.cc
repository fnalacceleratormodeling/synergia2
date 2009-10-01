#include "mpi_comm_wrap.h"
#include <stdexcept>

MPI_comm_wrap::MPI_comm_wrap(MPI_Comm comm)
{
    this->comm = comm;
}

int
MPI_comm_wrap::get_rank()
{
    int error, rank;
    error = MPI_Comm_rank(comm, &rank);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_rank");
    }
    return rank;
}

int
MPI_comm_wrap::get_size()
{
    int error, size;
    error = MPI_Comm_size(comm, &size);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_size");
    }
    return size;
}

MPI_Comm
MPI_comm_wrap::get()
{
    return comm;
}
