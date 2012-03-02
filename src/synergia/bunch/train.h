#ifndef TRAIN_H_
#define TRAIN_H_
#include <vector>

#include "synergia/bunch/bunch.h"

class Train_comms
{
private:
    int num_bunches;
    std::vector<bool > on_this_rank;
    Commxx master_comm;
    std::vector<Commxx > comms;
    MPI_Group master_group;
    std::vector<MPI_Group > groups;
    std::vector< int> proc_counts;
    std::vector< int> proc_offsets;
protected:
    void
    verify_index(int index) const;
public:
    /// @param num_bunches the number of bunches in the train
    /// @param bunch_separation the spatial separation between adjacent bunches [m]
    /// @param master_comm the communicator for the entire train
    Train_comms(int num_bunches, Commxx const& master_comm);
    /// Get the number of bunches in the train
    int
    get_num_bunches() const;
    /// Get the master communicator
    Commxx const&
    get_master_comm() const;
    /// Get the communicator for a bunch
    /// @param index which bunch
    Commxx const &
    get_comm(int index) const;
    /// Return true if bunch with index index overlaps with the current processor
    /// @param index which bunch
    bool
    is_on_this_rank(int index) const;

    std::vector< int> get_proc_counts() const;
    std::vector< int> get_proc_offsets() const;
    virtual
    ~Train_comms();
};

typedef boost::shared_ptr<Train_comms > Train_comms_sptr;

/// Represents a train of bunches separated by a fixed distance
class Bunch_train: public Train_comms
{
private:
   double bunch_separation;
   std::vector<Bunch_sptr > bunches;
public:
   /// @param num_bunches the number of bunches in the train
   /// @param bunch_separation the spatial separation between adjacent bunches [m]
   /// @param master_comm the communicator for the entire train
       Bunch_train(int num_bunches, double bunch_separation, Commxx const& master_comm);
   /// Get the spatial separation between adjacent bunches [m]
   double
   get_bunch_separation() const;
   /// Set the bunch with index index
   /// @param index which bunch
   void
   set_bunch_sptr(int index, Bunch_sptr bunch_sptr);
   /// Get the bunch with index index
   /// @param index which bunch
   Bunch_sptr
   get_bunch_sptr(int index);


};

typedef boost::shared_ptr<Bunch_train > Bunch_train_sptr;


#if 0
/// Represents a train of bunch_with_diagnostics separated by a fixed distance
class Bunch_with_diagnostics_train: public Train_comms
{
private:
    double bunch_separation;
    std::vector<Bunch_with_diagnostics_sptr > bunch_diags;
public:
        Bunch_with_diagnostics_train(int num_bunches, double bunch_separation, Commxx const& master_comm);
     /// Get the spatial separation between adjacent bunches [m]
    double
    get_bunch_separation() const;
    /// Set the bunch_diag with index index
    /// @param index which bunch_diag
    void
    set_bunch_diag_sptr(int index, Bunch_with_diagnostics_sptr bunch_diag_sptr);
    /// Get the bunch_diag with index index
    /// @param index which bunch
    Bunch_with_diagnostics_sptr
    get_bunch_diag_sptr(int index);

};

typedef boost::shared_ptr<Bunch_with_diagnostics_train > Bunch_with_diagnostics_train_sptr;
#endif

#endif /* BUNCH_DIAG_TRAIN_H_ */
