#include "synergia/bunch/train.h"
#include <stdexcept>
#include "synergia/utils/parallel_utils.h"

Train_comms::Train_comms(int num_bunches, const Commxx & master_comm) :
    num_bunches(num_bunches), on_this_rank(num_bunches), master_comm(master_comm),
    comms(num_bunches), groups(num_bunches),proc_counts(master_comm.get_size(),0),
    proc_offsets(master_comm.get_size(),0)
{
   

    int error;
    error = MPI_Comm_group(master_comm.get(), &master_group);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in Train(MPI_Comm_group(1))");
    }
    std::vector<std::vector<int > > ranks(
            distribute_1d(master_comm, num_bunches));
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        on_this_rank[bunch] = false;
        for (int i = 0; i < ranks[bunch].size(); ++i) {
            if (ranks[bunch][i] == master_comm.get_rank()) {
                on_this_rank[bunch] = true;
            }
        }
        error = MPI_Group_incl(master_group, ranks[bunch].size(),
                &ranks[bunch][0], &groups[bunch]);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Train(MPI_Group_incl)");
        }
        MPI_Comm bunch_comm;
        error = MPI_Comm_create(master_comm.get(), groups[bunch], &bunch_comm);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Train(MPI_Comm_create)");
        }
        comms[bunch].set(bunch_comm);
         
    }

     counts_and_offsets_for_impedance(master_comm,num_bunches,proc_offsets,proc_counts); 

}



int
Train_comms::get_num_bunches() const
{
    return num_bunches;
}


const Commxx &
Train_comms::get_master_comm() const
{
    return master_comm;
}

void
Train_comms::verify_index(int index) const
{
    if ((index < 0) || (index >= num_bunches)) {
        throw std::runtime_error("Bunch_with_diagnostics_train: invalid index");
    }
}

const Commxx &
Train_comms::get_comm(int index) const
{
    verify_index(index);
    return comms[index];
}

    
bool
Train_comms::is_on_this_rank(int index) const
{
    verify_index(index);
    return on_this_rank[index];
}


std::vector< int> 
Train_comms::get_proc_counts() const
{
 return proc_counts;
}

std::vector< int> 
Train_comms::get_proc_offsets() const
{
return proc_offsets;
}


Train_comms::~Train_comms()
{
    int error;
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        if (on_this_rank[bunch]) {
            MPI_Comm comm(comms[bunch].get());
            error = MPI_Comm_free(&comm);
            if (error != MPI_SUCCESS) {
                throw std::runtime_error(
                        "MPI error in Bunch_train(MPI_Comm_free)");
            }
        }
        error = MPI_Group_free(&groups[bunch]);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in Bunch_train(MPI_Group_free)");
        }
    }
}


//*******************************************************************************

Bunch_train::Bunch_train(int num_bunches, double bunch_separation,
        const Commxx & master_comm) : Train_comms(num_bunches, master_comm), 
        bunches(num_bunches), bunch_separation(bunch_separation)
{
}        

double
Bunch_train::get_bunch_separation() const
{
    return bunch_separation;
}

void
Bunch_train::set_bunch_sptr(int index, Bunch_sptr bunch_sptr)
{
    verify_index(index);
    bunches[index] = bunch_sptr;
}

Bunch_sptr
Bunch_train::get_bunch_sptr(int index)
{
    verify_index(index);
     return bunches[index];
}


//*******************************************************************************


Bunch_with_diagnostics_train::Bunch_with_diagnostics_train(int num_bunches, double bunch_separation,
        const Commxx & master_comm) : Train_comms(num_bunches,  master_comm),
        bunch_diags(num_bunches), bunch_separation(bunch_separation)
{
}        


double
Bunch_with_diagnostics_train::get_bunch_separation() const
{
    return bunch_separation;
}
void
Bunch_with_diagnostics_train::set_bunch_diag_sptr(int index, Bunch_with_diagnostics_sptr bunch_diag_sptr)
{
    verify_index(index);
    bunch_diags[index] = bunch_diag_sptr;
}

Bunch_with_diagnostics_sptr
Bunch_with_diagnostics_train::get_bunch_diag_sptr(int index)
{
     verify_index(index);
     return bunch_diags[index];
}
