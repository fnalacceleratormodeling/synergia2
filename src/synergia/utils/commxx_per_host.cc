#include "commxx_per_host.h"

#include <limits.h>
#include <stdexcept>

static size_t hash(const char * s)
{
  return (*(size_t*)s) >> 2;
}

Commxx_per_host::Commxx_per_host()
: Commxx()
{
  char name[MPI_MAX_PROCESSOR_NAME];
  int  name_len;
  MPI_Get_processor_name( name, &name_len );

  int color = hash(name) % INT_MAX;

  int result = MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);
  if( result!=MPI_SUCCESS ) throw std::runtime_error("MPI error in MPI_Comm_split");
}
