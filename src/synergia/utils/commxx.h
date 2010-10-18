#ifndef COMMXX_H_
#define COMMXX_H_

#include "mpi.h"

/// Commxx is a thin wrapper around MPI communicator (MPI_Comm) objects.
/// This is a C++-only class.
///
/// Python: The equivalent functionality is provided by mpi4py Comm objects.
/// mpi4py Comm objects may be passed from python anywhere a Commxx object is
/// expected.
class Commxx
{
private:
    MPI_Comm comm;
public:
    /// Construct a Commxx object from an MPI_Comm communicator
    /// @param comm is an MPI communicator, e.g., MPI_COMM_WORLD
    Commxx(MPI_Comm comm);

    /// Construct a Commxx object using MPI_COMM_WORLD
    Commxx();

    /// Get communicator rank.
    int
    get_rank() const;

    /// Get communicator size
    int
    get_size() const;

    /// Change the underlying communicator
    /// @param comm is an MPI communicator, e.g., MPI_COMM_WORLD
    void
    set(MPI_Comm comm);

    /// Extract the MPI_comm object wrapped by the Commxx instance.
    MPI_Comm
    get() const;
};

#endif /* COMMXX_H_ */
