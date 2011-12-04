#ifndef COMMXX_H_
#define COMMXX_H_

#include <stdexcept>
#include "mpi.h"

#include "synergia/utils/serialization.h"

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

    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const
        {
            int result;
            MPI_Comm_compare(MPI_COMM_WORLD, comm, &result);
            if (result != MPI_IDENT) {
                throw std::runtime_error(
                        "Commxx: can only serialize MPI_COMM_WORLD");
            }
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
        {
            comm = MPI_COMM_WORLD;
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

#endif /* COMMXX_H_ */
