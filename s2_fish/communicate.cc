#include "communicate.h"
#include <mpi.h>

void
gather_rho(Real_scalar_field &rho, int upper_limit)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int uppers[size];
    int lengths[size];
    MPI_Allgather(reinterpret_cast<void*>(&upper_limit), 1, MPI_INT,
                  reinterpret_cast<void*>(uppers), 1, MPI_INT, MPI_COMM_WORLD);
    Int3 shape(rho.get_points().get_shape());
    lengths[0] = uppers[0] * shape[1] * shape[2];
    for (int i = 1; i < size; ++i) {
        if (uppers[i] > shape[0]) {
            uppers[i] = shape[0];
        }
        lengths[i] = (uppers[i] - uppers[i-1]) * shape[1] * shape[2];
        if (lengths[i] < 0) {
            lengths[i] = 0;
        }
    }
    void *tmp;
    if (lengths[rank] > 0) {
        tmp = malloc(lengths[rank] * sizeof(double));
    }
    MPI_Reduce_scatter(reinterpret_cast<void*>
                       (rho.get_points().get_base_address()),
                       tmp, lengths, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int offset;
    if (rank == 0) {
        offset = 0;
    } else {
        offset = uppers[rank-1];
    }
    if ((offset < shape[0]) && (lengths[rank] > 0)) {
        memcpy(reinterpret_cast<void*>
               (rho.get_points().get_offset_base_address(offset)),
               tmp, lengths[rank]*sizeof(double));
    }
    if (lengths[rank] > 0) {
        free(tmp);
    }
}

void
gather_global_rho(Real_scalar_field &local_rho,Real_scalar_field &global_rho)
{
    MPI_Allreduce(reinterpret_cast<void*>(local_rho.get_points().get_base_address()),
                reinterpret_cast<void*>(global_rho.get_points().get_base_address()),
                local_rho.get_points().get_length(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void
allreduce_broadcast_E(Real_scalar_field &E, int i_lower, int i_upper)
{

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size_t byte_length = E.get_points().get_length() * sizeof(double);
    void *buffer = malloc(byte_length);
    MPI_Allreduce(reinterpret_cast<void*>(E.get_points().get_base_address()),
                  buffer, E.get_points().get_length(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(reinterpret_cast<void*>(E.get_points().get_base_address()),
           buffer, byte_length);
    int uppers[size];
    int lengths[size];
    int displacements[size];
    MPI_Allgather(reinterpret_cast<void*>(&i_upper), 1, MPI_INT,
                  reinterpret_cast<void*>(uppers), 1, MPI_INT, MPI_COMM_WORLD);
    free(buffer);
}

void
old_broadcast_E(Real_scalar_field &E, int i_lower, int i_upper)
{
    //~ int rank,size;
    //~ MPI_Comm_size(MPI_COMM_WORLD, &size);
    //~ if (size == 0) {
    //~ return;
    //~ }
    //~ MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //~ int uppers[size];
    //~ MPI_Allgather(reinterpret_cast<void*>(&i_upper),1,MPI_INT,
    //~ reinterpret_cast<void*>(uppers),1,MPI_INT,MPI_COMM_WORLD);
    //~ Int3 shape(E.get_points().get_shape());
    //~ int last_upper = 0;
    //~ int sender = 0;
    //~ while (last_upper<shape[0]) {
    //~ MPI_Bcast(reinterpret_cast<void*>(E.get_points().get_offset_base_address(last_upper)),
    //~ (uppers[sender]-last_upper)*shape[1]*shape[2],
    //~ MPI_DOUBLE, sender, MPI_COMM_WORLD);
    //~ last_upper = uppers[sender];
    //~ sender++;
    //~ }

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int uppers[size];
    int lengths[size];
    int displacements[size];
    MPI_Allgather(reinterpret_cast<void*>(&i_upper), 1, MPI_INT,
                  reinterpret_cast<void*>(uppers), 1, MPI_INT, MPI_COMM_WORLD);
    Int3 shape(E.get_points().get_shape());
    lengths[0] = uppers[0] * shape[1] * shape[2];
    displacements[0] = 0;
    for (int i = 1; i < size; ++i) {
        displacements[i] = uppers[i-1] * shape[1] * shape[2];
        if (uppers[i] > shape[0]) {
            uppers[i] = shape[0];
        }
        lengths[i] = (uppers[i] - uppers[i-1]) * shape[1] * shape[2];
    }
    int offset;
    if (rank == 0) {
        offset = 0;
    } else {
        offset = uppers[rank-1];
    }
    void* tmp = malloc(E.get_points().get_length() * sizeof(double));
    MPI_Allgatherv(reinterpret_cast<void*>
                   (E.get_points().get_offset_base_address(offset)),
                   lengths[rank], MPI_DOUBLE,
                   tmp, lengths, displacements, MPI_DOUBLE, MPI_COMM_WORLD);
    memcpy(reinterpret_cast<void*>
           (E.get_points().get_base_address()),
           tmp, E.get_points().get_length()*sizeof(double));
    free(tmp);
}

void
allgather_phi(Real_scalar_field &local_phi, Real_scalar_field &full_phi)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        full_phi.get_points().copy(&local_phi.get_points());
        return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int upper_limit = local_phi.get_points().get_dim0_upper();
    int uppers[size];
    int lowers[size];
    int lengths[size];
    int num_with_data = 1;
    MPI_Allgather(reinterpret_cast<void*>(&upper_limit), 1, MPI_INT,
                  reinterpret_cast<void*>(uppers), 1, MPI_INT, MPI_COMM_WORLD);
    Int3 shape(full_phi.get_points().get_shape());
    lengths[0] = uppers[0] * shape[1] * shape[2];
    lowers[0] = 0;
    for (int i = 1; i < size; ++i) {
        lowers[i] = uppers[i-1];
        if (lowers[i] < shape[0]) {
            ++num_with_data;
        }
        if (uppers[i] > shape[0]) {
            uppers[i] = shape[0];
        }
        lengths[i] = (uppers[i] - uppers[i-1]) * shape[1] * shape[2];
    }
    MPI_Group group_world, group_half;
    MPI_Comm comm_half;
    int ranks[num_with_data];
    for (int i = 0; i < num_with_data; ++i) {
        ranks[i] = i;
    }
    //~ if (rank == 0) {
    //~ std::cout << "jfa: num_with_data = " << num_with_data << std::endl;
    //~ int total=0;
    //~ for(int i=0; i<num_with_data; ++i) {
    //~ std::cout << i << ": " << lowers[i] << " " << lengths[i] << std::endl;
    //~ total += lengths[i];
    //~ }
    //~ std::cout << "total lengths = " << total << ", phi length = " << full_phi.get_points().get_length() << std::endl;
    //~ }
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, num_with_data, ranks, &group_half);
    MPI_Comm_create(MPI_COMM_WORLD, group_half, &comm_half);
    if (rank < num_with_data) {
        MPI_Gatherv(reinterpret_cast<void*>
                    (local_phi.get_points().get_offset_base_address(lowers[rank])),
                    lengths[rank], MPI_DOUBLE,
                    reinterpret_cast<void*>
                    (full_phi.get_points().get_base_address()), lengths, lowers,
                    MPI_DOUBLE, 0, comm_half);
    }
    MPI_Bcast(reinterpret_cast<void*>
              (full_phi.get_points().get_base_address()),
              full_phi.get_points().get_length(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //~ MPI_Comm_free(&comm_half);
}

void
fill_guards(Real_scalar_field &phi, Fftw_helper &fftwh)
{
    Int3 shape(phi.get_points().get_shape());
    size_t message_size = shape[1] * shape[2];
    void *recv_buffer, *send_buffer;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    if (size < 3) {
        return;
    }
    if (fftwh.lower() >= phi.get_points().get_shape()[0]) {
        return;
    }
    // send to the right
    recv_buffer = reinterpret_cast<void*>(phi.get_points().get_offset_base_address(fftwh.guard_lower()));
    send_buffer = reinterpret_cast<void*>(phi.get_points().get_offset_base_address(fftwh.upper() - 1));
    if (fftwh.upper() != fftwh.guard_upper()) {
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
    }
    if (fftwh.lower() != fftwh.guard_lower()) {
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
    }

    //send to the left
    recv_buffer = reinterpret_cast<void*>(phi.get_points().get_offset_base_address(fftwh.guard_upper() - 1));
    send_buffer = reinterpret_cast<void*>(phi.get_points().get_offset_base_address(fftwh.lower()));
    if (fftwh.lower() != fftwh.guard_lower()) {
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
    }
    if (fftwh.upper() != fftwh.guard_upper()) {
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
    }
}

void
broadcast_E(Real_scalar_field &E, int i_lower, int i_upper)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        return;
    }
    int uppers[size];
    int lowers[size];
    int lengths[size];
    int num_with_data = 1;
    MPI_Allgather(reinterpret_cast<void*>(&i_upper), 1, MPI_INT,
                  reinterpret_cast<void*>(uppers), 1, MPI_INT, MPI_COMM_WORLD);
    Int3 shape(E.get_points().get_shape());
    lengths[0] = uppers[0] * shape[1] * shape[2];
    lowers[0] = 0;
    for (int i = 1; i < size; ++i) {
        lowers[i] = uppers[i-1] * shape[1] * shape[2];
        if (lowers[i] < shape[0]*shape[1]*shape[2]) {
            ++num_with_data;
        }
        if (uppers[i] > shape[0]) {
            uppers[i] = shape[0];
        }
        lengths[i] = (uppers[i] - uppers[i-1]) * shape[1] * shape[2];
    }
    MPI_Group group_world, group_half;
    MPI_Comm comm_half;
    int ranks[num_with_data];
    for (int i = 0; i < num_with_data; ++i) {
        ranks[i] = i;
    }
    //~ if (rank == 0) {
    //~ std::cout << "jfa: num_with_data = " << num_with_data << std::endl;
    //~ int total=0;
    //~ for(int i=0; i<num_with_data; ++i) {
    //~ std::cout << i << ": " << lowers[i] << " " << lengths[i] << std::endl;
    //~ total += lengths[i];
    //~ }
    //~ std::cout << "total lengths = " << total << ", phi length = " << E.get_points().get_length() << std::endl;
    //~ }
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, num_with_data, ranks, &group_half);
    MPI_Comm_create(MPI_COMM_WORLD, group_half, &comm_half);
    void *buffer = malloc(E.get_points().get_length() * sizeof(double));
    if (rank < num_with_data) {
        MPI_Gatherv(reinterpret_cast<void*>
                    (E.get_points().get_base_address() + lowers[rank]),
                    lengths[rank], MPI_DOUBLE,
                    buffer, lengths, lowers,
                    MPI_DOUBLE, 0, comm_half);
    }
    MPI_Bcast(buffer,
              E.get_points().get_length(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    memcpy(reinterpret_cast<void*>(E.get_points().get_base_address()),
           buffer, E.get_points().get_length()*sizeof(double));
    free(buffer);
    //~ MPI_Comm_free(&comm_half);
}
// Collect the stripes, and fill the potential over the entire physical size 
// of the problem, on every node. 
void
broadcast_Phi(Real_scalar_field &E, Real_scalar_field &EAll, int i_lower, int i_upper)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        memcpy(reinterpret_cast<void*>(EAll.get_points().get_base_address()),
	       reinterpret_cast<void*>(E.get_points().get_base_address()),
	       E.get_points().get_length()*sizeof(double));
        return;
    }
    int uppers[size];
    int lowersSend[size];  // on the sending node 
    int lowersRcv[size];  // on the receiving node (rank zero, or root) 
    int lengths[size];
    int num_with_data = 1;
    Int3 shape(E.get_points().get_shape());
    int i_upperEff = i_upper-1; // Skip the guards 
    if (i_upper == shape[0]) i_upperEff=i_upper; // except at the upper edge of entire grid
    MPI_Allgather(reinterpret_cast<void*>(&i_upperEff), 1, MPI_INT,
                  reinterpret_cast<void*>(uppers), 1, MPI_INT, MPI_COMM_WORLD);
    lengths[0] = uppers[0] * shape[1] * shape[2];
    lowersRcv[0] = 0;
    lowersSend[0] = 0;
    for (int i = 1; i < size; ++i) {
        lowersRcv[i] = uppers[i-1] * shape[1] * shape[2];
        lowersSend[i] = shape[1] * shape[2];
        if (lowersRcv[i] < shape[0]*shape[1]*shape[2]) {
            ++num_with_data;
        }
        if (uppers[i] > shape[0]) {
            uppers[i] = shape[0];
        }
        lengths[i] = (uppers[i] - uppers[i-1]) * shape[1] * shape[2];
    }
//    std::cerr << " ... Upper/lower recalc, rank " << rank << " --- ";
//    for (int kS=0; kS != size; kS++) 
//      std::cerr << " r- " << kS << " lowR=" << lowersRcv[kS] << " lowS=" << 
//                << " upp=" 
//                << uppers[kS] << " len=" << lengths[kS];
//     std::cerr << std::endl;
    MPI_Group group_world, group_half;
    MPI_Comm comm_half;
    int ranks[num_with_data];
    for (int i = 0; i < num_with_data; ++i) {
        ranks[i] = i;
    }
    //~ if (rank == 0) {
    //~ std::cout << "jfa: num_with_data = " << num_with_data << std::endl;
    //~ int total=0;
    //~ for(int i=0; i<num_with_data; ++i) {
    //~ std::cout << i << ": " << lowers[i] << " " << lengths[i] << std::endl;
    //~ total += lengths[i];
    //~ }
    //~ std::cout << "total lengths = " << total << ", phi length = " << E.get_points().get_length() << std::endl;
    //~ }
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, num_with_data, ranks, &group_half);
    MPI_Comm_create(MPI_COMM_WORLD, group_half, &comm_half);
    void *buffer = malloc(EAll.get_points().get_length() * sizeof(double));
    if (rank < num_with_data) {
        
//        double *aValTmp=E.get_points().get_base_address();
//	if (rank == 1) { 
//            aValTmp=E.get_points().get_base_address();
//	    for (int k=0; k != E.get_points().get_length(); k++) 
//	      std::cerr << " Before GatherV, from rank " << rank << 
//	     " .. k = " << k << " val = " << aValTmp[k] << std::endl;
//	}
        MPI_Gatherv(reinterpret_cast<void*>
                    (E.get_points().get_base_address() + lowersSend[rank]),
                    lengths[rank], MPI_DOUBLE,
                    buffer, lengths, lowersRcv,
                    MPI_DOUBLE, 0, comm_half);
	  // We only have values for rank 0... (root) 
//	  if (rank == 0) { 
//            aValTmp=reinterpret_cast<double*> (buffer);
//	    for (int k=0; k != EAll.get_points().get_length(); k++) 
//	      std::cerr << " After GatherV, from rank " << rank << 
//	     " .. k = " << k << " val = " << aValTmp[k] << std::endl;
//	  }
		    
    }
    
    MPI_Bcast(buffer,
              EAll.get_points().get_length(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    memcpy(reinterpret_cast<void*>(EAll.get_points().get_base_address()),
           buffer, EAll.get_points().get_length()*sizeof(double));
    free(buffer);
    //~ MPI_Comm_free(&comm_half);
}
