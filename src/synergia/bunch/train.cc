#include "synergia/bunch/train.h"
#include <stdexcept>
#include "synergia/utils/parallel_utils.h"

Train_comms::Train_comms(int num_bunches, Commxx_sptr master_comm_sptr) :
    num_bunches(num_bunches), master_comm_sptr(master_comm_sptr),
    comms(0), proc_counts(master_comm_sptr->get_size(),0),
    proc_offsets(master_comm_sptr->get_size(),0)
{


    std::vector<std::vector<int > > ranks(
            distribute_1d(*master_comm_sptr, num_bunches));
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        comms.push_back(Commxx_sptr(new Commxx(master_comm_sptr,ranks[bunch])));
    }

     counts_and_offsets_for_impedance(*master_comm_sptr,num_bunches,proc_offsets,proc_counts);

}



int
Train_comms::get_num_bunches() const
{
    return num_bunches;
}


const Commxx &
Train_comms::get_master_comm() const
{
    return *master_comm_sptr;
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
    return *comms[index];
}

Commxx_sptr
Train_comms::get_comm_sptr(int index) const
{
    verify_index(index);
    return comms[index];
}


bool
Train_comms::is_on_this_rank(int index) const
{
    verify_index(index);
    return comms[index]->has_this_rank();
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
}


//*******************************************************************************

Bunch_train::Bunch_train(int num_bunches, double bunch_separation,
        Commxx_sptr master_comm_sptr) : Train_comms(num_bunches, master_comm_sptr),
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
    if (bunch_sptr->get_comm().get() != get_comm(index).get()) throw
        std::runtime_error("Bunch_train::set_bunch_sptr-- the train and the bunch communicators are different");
    bunches[index] = bunch_sptr;
}

Bunch_sptr
Bunch_train::get_bunch_sptr(int index)
{
    verify_index(index);
     return bunches[index];
}


//*******************************************************************************

#if 0
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
    bunch_diag_sptr->check_bunch_pointer_in_diagnostics();
    if (bunch_diag_sptr->get_comm().get() != get_comm(index).get()) throw
        std::runtime_error("Bunch_with_diagnostics_train::set_bunch_diag_sptr-- the train and the bunch communicators are different");
    bunch_diags[index] = bunch_diag_sptr;
}

Bunch_with_diagnostics_sptr
Bunch_with_diagnostics_train::get_bunch_diag_sptr(int index)
{
     verify_index(index);
     return bunch_diags[index];
}
#endif
