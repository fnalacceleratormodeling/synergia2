#ifndef BUNCH_TRAIN_H_
#define BUNCH_TRAIN_H_

#include "synergia/bunch/bunch.h"

class Bunch_train
{
private:
    Bunches bunches;
    std::vector<double > spacings;      
    bool has_parent_comm;
    Commxx_sptr parent_comm_sptr;
    /// counts and offsets are needed for impedance, counts.size()=offsets.size()=num_procs
    /// they are meaningfull only on the local rank=0 of every bunch communicator 
    /// for example:  3 bunches on 5 processors: proc 2 and 4 has local rank different form zero
    ///              (count[0]=1, offset[0]=0), (count[1]=1, offset[3]=1), (count[3]=1, offset[3]=2)
    /// for example:  5 bunches on 3 processors: 
    ///                (count[0]=1, offset[0]=0), (count[1]=2, offset[3]=1), (count[2]=2, offset[2]=3)
    std::vector< int> proc_counts_imped; 
    std::vector< int> proc_offsets_imped; 
    void
    set_bucket_indices();
    void find_parent_comm_sptr();
    void calculates_counts_and_offsets_for_impedance();
public:
    Bunch_train(Bunches const& bunches, double spacing);
    Bunch_train(Bunches const& bunches, std::vector<double > const& spacings);
    // Default constructor for serialization use only
    Bunch_train();
    Commxx_sptr 
    get_parent_comm_sptr();
    size_t
    get_size() const;
    Bunches &
    get_bunches();
    std::vector<double > &
    get_spacings();
    
    // update the total particle number for all bunches in the bunch train
    // note that calling each bunch's update_total_num() wont do the actual
    // work if the caller's rank is not part of the bunch's communicator.
    // on the other hand, calling bunch_train's update_bunch_total_num() 
    // gurantees that all bunches are updated for all ranks
    void 
    update_bunch_total_num();
    
    std::vector< int> &
    get_proc_counts_for_impedance();
    std::vector< int> &
    get_proc_offsets_for_impedance();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Bunch_train();
};

typedef boost::shared_ptr<Bunch_train > Bunch_train_sptr;

#endif /* BUNCH_TRAIN_H_ */
