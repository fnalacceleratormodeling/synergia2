#include "bunch_train.h"
#include <stdexcept>
#include "synergia/utils/parallel_utils.h"




Bunch_train::Bunch_train(int num_bunches, double bunch_separation,
        const Commxx & master_comm) :
    num_bunches(num_bunches), bunch_separation(bunch_separation),
    bunches(num_bunches), on_this_rank(num_bunches), master_comm(master_comm),
    comms(num_bunches), groups(num_bunches)
{
   
    int error;
    error = MPI_Comm_group(master_comm.get(), &master_group);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in Bunch_train(MPI_Comm_group(1))");
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
                    "MPI error in Bunch_train(MPI_Group_incl)");
        }
        MPI_Comm bunch_comm;
        error = MPI_Comm_create(master_comm.get(), groups[bunch], &bunch_comm);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Bunch_train(MPI_Comm_create)");
        }
        comms[bunch].set(bunch_comm);
    }
}

int
Bunch_train::get_num_bunches() const
{
    return num_bunches;
}

double
Bunch_train::get_bunch_separation() const
{
    return bunch_separation;
}

const Commxx &
Bunch_train::get_master_comm() const
{
    return master_comm;
}

void
Bunch_train::verify_index(int index) const
{
    if ((index < 0) || (index >= num_bunches)) {
        throw std::runtime_error("Bunch_train: invalid index");
    }
}

const Commxx &
Bunch_train::get_comm(int index) const
{
    verify_index(index);
    return comms[index];
}

bool
Bunch_train::is_on_this_rank(int index) const
{
    verify_index(index);
    return on_this_rank[index];
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

Bunch_train::~Bunch_train()
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
