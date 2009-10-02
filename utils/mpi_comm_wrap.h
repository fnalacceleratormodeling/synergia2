#ifndef MPI_COMM_WRAP_H_
#define MPI_COMM_WRAP_H_

#include "mpi.h"

class MPI_comm_wrap
{
private:
    MPI_Comm comm;
public:
    MPI_comm_wrap(MPI_Comm comm);
    int
    get_rank() const;
    int
    get_size() const;
    MPI_Comm
    get();
};

#endif /* MPI_COMM_WRAP_H_ */
