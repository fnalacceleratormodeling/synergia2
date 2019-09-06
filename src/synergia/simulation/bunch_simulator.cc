#include <numeric>

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/independent_operation.h"


Bunch_simulator
Bunch_simulator::create_single_bunch_simulator(
        Reference_particle const & ref,
        size_t num_particles,
        double num_real_particles,
        Commxx const & comm )
{
    return construct( ref, ref, 
                      num_particles, num_real_particles,
                      1, 0,     // num_bunches
                      1.0, 1.0, // spacing
                      comm );
}


Bunch_simulator 
Bunch_simulator::create_bunch_train_simulator(
        Reference_particle const & ref,
        size_t num_particles,
        double num_real_particles,
        size_t num_bunches,
        double spacing,
        Commxx const & comm )
{
    return construct( ref, ref, 
                      num_particles, num_real_particles,
                      num_bunches, 0, // num_bunches
                      spacing, 1.0,   // spacing
                      comm );
}

Bunch_simulator
Bunch_simulator::construct(
        Reference_particle const & ref_pri,
        Reference_particle const & ref_sec,
        size_t num_part,
        double num_real_part,
        size_t num_bunches_pri,
        size_t num_bunches_sec,
        double spacing_pri,
        double spacing_sec,
        Commxx const & comm)
{
    int size = comm.size();
    int rank = comm.rank();

    size_t num_bunches = num_bunches_pri + num_bunches_sec;

    Commxx comm_pri;
    Commxx comm_sec;

    std::vector<std::vector<int>> bunch_ranks;

    if ( size < num_bunches )
    {
        if (num_bunches % size != 0)
        {
            throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of bunches must be divisible by the number of ranks" );
        }

        size_t bunches_per_rank = num_bunches / size;

        size_t num_ranks_pri = std::ceil(1.0*num_bunches_pri/bunches_per_rank);
        size_t num_ranks_sec = std::ceil(1.0*num_bunches_sec/bunches_per_rank);

        size_t st_starting_rank = size - num_ranks_sec;

        // which ranks (in the context of simulator comm) 
        // are a particular bunch at
        for (int b=0; b<num_bunches; ++b)
        {
            int rank = b / bunches_per_rank;
            bunch_ranks.emplace_back(1, rank);
        }

        // construct the communicators for the 
        // primary and secondary trains
        std::vector<int> p_ranks(num_ranks_pri);
        std::vector<int> s_ranks(num_ranks_sec);

        for (int r=0; r<num_ranks_pri; ++r) p_ranks[r] = r;
        comm_pri = comm.group(p_ranks);

        for (int r=0; r<num_ranks_sec; ++r) s_ranks[r] = st_starting_rank + r;
        comm_sec = comm.group(s_ranks);
    }
    else
    {
        if (size % num_bunches != 0)
        {
            throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of bunches must be divisible by the number of ranks" );
        }

        size_t ranks_per_bunch = size / num_bunches;

        size_t num_ranks_pri = ranks_per_bunch * num_bunches_pri;
        size_t num_ranks_sec = ranks_per_bunch * num_bunches_sec;

        size_t st_starting_rank = num_ranks_pri;

        // which ranks (in the context of simulator comm) 
        // do a particular bunch resides
        for (int b=0; b<num_bunches; ++b)
        {
            int starting_rank = b * ranks_per_bunch;
            std::vector<int> ranks(ranks_per_bunch);
            std::iota(ranks.begin(), ranks.end(), starting_rank);
            bunch_ranks.emplace_back(std::move(ranks));
        }

        // construct the communicators for the 
        // primary and secondary trains
        std::vector<int> p_ranks(num_ranks_pri);
        std::vector<int> s_ranks(num_ranks_sec);

        for (int r=0; r<num_ranks_pri; ++r) p_ranks[r] = r;
        comm_pri = comm.group(p_ranks);

        for (int r=0; r<num_ranks_sec; ++r) s_ranks[r] = st_starting_rank + r;
        comm_sec = comm.group(s_ranks);
    }

    return Bunch_simulator( Bunch_train( ref_pri, 
                                         num_bunches_pri, 
                                         num_part, 
                                         num_real_part, 
                                         spacing_pri, 
                                         comm_pri ),
                            Bunch_train( ref_sec, 
                                         num_bunches_sec, 
                                         num_part, 
                                         num_real_part, 
                                         spacing_sec, 
                                         comm_sec ),
                            bunch_ranks, 
                            comm );
}




Bunch_simulator::Bunch_simulator(
        Bunch_train && pt, 
        Bunch_train && st,
        std::vector<std::vector<int>> const & bunch_ranks,
        Commxx const & comm )
    : trains{std::move(pt), std::move(st)}
    , bunch_ranks(bunch_ranks)
    , pt_bunches(pt.get_size())
    , st_bunches(st.get_size())
    , comm(comm)
    , diags_step()
    , diags_loss()
    , diags_ele()
    , diags_opr()
    , diags_opn()
{
}

void
Bunch_simulator::diag_action_step_and_turn(int turn_num, int step_num)
{
    for (auto const & dt : diags_step)
    {
        if (dt.trigger(turn_num, step_num))
        {
            get_bunch(dt.train, dt.bunch)
                .diag_update_and_write(dt.diag_name);
        }
    }
}

void
Bunch_simulator::diag_action_particle_loss_update()
{
}

void
Bunch_simulator::diag_action_particle_loss_write()
{
}

void
Bunch_simulator::diag_action_element(Lattice_element const & element)
{
}


void
Bunch_simulator::diag_action_operator(Operator const & opr)
{
}

void
Bunch_simulator::diag_action_operation(Independent_operation const & opn)
{
}

void
Bunch_simulator::prop_action_first(Lattice & lattice)
{
}

void
Bunch_simulator::prop_action_step_end(Lattice & lattice, int turn, int step)
{
}

void
Bunch_simulator::prop_action_turn_end(Lattice & lattice, int turn)
{
}

void 
Bunch_simulator::set_lattice_reference_particle(Reference_particle const & ref)
{
}
