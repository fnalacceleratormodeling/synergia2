#ifndef COMMXX_H_
#define COMMXX_H_

#include "mpi.h"

class Commxx
{
private:
    MPI_Comm comm;
public:
    Commxx(MPI_Comm comm);
    int
    get_rank() const;
    int
    get_size() const;
    MPI_Comm
    get() const;
};

#endif /* COMMXX_H_ */
