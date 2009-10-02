#include "commxx.h"
#include <stdexcept>

Commxx::Commxx(MPI_Comm comm)
{
    this->comm = comm;
}

int
Commxx::get_rank() const
{
    int error, rank;
    error = MPI_Comm_rank(comm, &rank);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_rank");
    }
    return rank;
}

int
Commxx::get_size() const
{
    int error, size;
    error = MPI_Comm_size(comm, &size);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_size");
    }
    return size;
}

MPI_Comm
Commxx::get() const
{
    return comm;
}
