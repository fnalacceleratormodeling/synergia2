#ifndef BUNCH_TRAIN_H_
#define BUNCH_TRAIN_H_
#include "synergia/bunch/bunch.h"
#include <vector>


/// Represents a train of bunches separated by a fixed distance
class Bunch_train
{
private:
    int num_bunches;
    double bunch_separation;
    std::vector<Bunch_sptr > bunches;
    std::vector<bool > on_this_rank;
    Commxx master_comm;
    std::vector<Commxx > comms;
    MPI_Group master_group;
    std::vector<MPI_Group > groups;
    void
    verify_index(int index) const;
public:
    /// @param num_bunches the number of bunches in the train
    /// @param bunch_separation the spatial separation between adjacent bunches [m]
    /// @param master_comm the communicator for the entire train    
        Bunch_train(int num_bunches, double bunch_separation, Commxx const& master_comm);
    /// Get the number of bunches in the train
    int
    get_num_bunches() const;
    /// Get the spatial separation between adjacent bunches [m]
    double
    get_bunch_separation() const;
    /// Get the master communicator
    Commxx const&
    get_master_comm() const;
    /// Get the communicator for a bunch
    /// @param index which bunch
    Commxx const&
    get_comm(int index) const;
    /// Return true if bunch with index index overlaps with the current processor
    /// @param index which bunch
    bool
    is_on_this_rank(int index) const;
    /// Set the bunch with index index
    /// @param index which bunch
    void
    set_bunch_sptr(int index, Bunch_sptr bunch_sptr);
    /// Get the bunch with index index
    /// @param index which bunch
    Bunch_sptr
    get_bunch_sptr(int index);
    
    

    // n.b. some communication routines left to be added

    ~Bunch_train();
};

typedef boost::shared_ptr<Bunch_train > Bunch_train_sptr;

#endif /* BUNCH_TRAIN_H_ */
