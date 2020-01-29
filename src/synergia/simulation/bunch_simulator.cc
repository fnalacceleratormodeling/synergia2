#include <numeric>

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/independent_operation.h"


Bunch_simulator
Bunch_simulator::create_empty_bunch_simulator()
{
    return construct( Reference_particle(), 
                      Reference_particle(),
                      1.0, 1.0,  // num_real_particle
                      0,   0,    // num bunches
                      1.0, 1.0,  // spacing
                      Commxx() );
}

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

    int num_ranks_pri = 0;
    int num_ranks_sec = 0;

    // in the case of num_bunches=0, both num_ranks_pri/sec would be 0, so
    // comm_pri/sec are both null communicators 
    if ( size < num_bunches )
    {
        if (num_bunches % size != 0)
        {
            throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of bunches must be divisible by the number of ranks" );
        }

        int bunch_per_rank = num_bunches / size;

        if (num_bunches_pri % bunch_per_rank != 0
                || num_bunches_sec % bunch_per_rank != 0 )
        {
            throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of bunches in primary or secondary train must be "
                    "divisible by the number of bunches per rank" );
        }

        num_ranks_pri = num_bunches_pri/bunch_per_rank;
        num_ranks_sec = num_bunches_sec/bunch_per_rank;
    }
    else if ( num_bunches > 0)
    {
        if (size % num_bunches != 0)
        {
            throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of ranks must be divisible by the number of bunches" );
        }

        int rank_per_bunch = size / num_bunches;

        num_ranks_pri = rank_per_bunch * num_bunches_pri;
        num_ranks_sec = rank_per_bunch * num_bunches_sec; 
    }

    // construct the communicators for the 
    // primary and secondary trains
    std::vector<int> p_ranks(num_ranks_pri);
    std::vector<int> s_ranks(num_ranks_sec);

    for (int r=0; r<num_ranks_pri; ++r) p_ranks[r] = r;
    comm_pri = comm.group(p_ranks);

    for (int r=0; r<num_ranks_sec; ++r) s_ranks[r] = num_ranks_pri + r;
    comm_sec = comm.group(s_ranks);

    return Bunch_simulator( Bunch_train( ref_pri, 
                                         num_bunches_pri, 
                                         num_part, 
                                         num_real_part, 
                                         spacing_pri, 
                                         comm_pri,
                                         0 ),
                            Bunch_train( ref_sec, 
                                         num_bunches_sec, 
                                         num_part, 
                                         num_real_part, 
                                         spacing_sec, 
                                         comm_sec,
                                         1 ),
                            comm );



    // A more general approach to the above logic could be like:
    //
    // bunch_ranks[ti][bi][ri]:
    //
    //   stores the MPI rank in the simulator communicator of given
    //   bunch index and rank inex
    //
    //   ti - index of the train (0 or 1)
    //   bi - index of the bunch (0 to num_bunches of the train)
    //   ri - index from 0 to number of ranks of this bunch
    //
    // bunch_idx_map[ti][bi]
    //
    //   stores the array index in the local bunch array of the
    //   given train/bunch. -1 if the bunch is not present in the
    //   current rank
    //
    //   ti: index of the train, 0 or 1
    //   bi: index of the bunch, 0 to num_bunches in the train
    //
    // E.g.: 
    //
    //   1. num_bunches_pri = 1, num_bunches_sec = 1
    //      total 2 bunches, 8 ranks (n_ranks > n_bunches)
    //
    //      bunch_ranks = { 
    //          { {0, 1, 2, 3} }, 
    //          { {4, 5, 6, 7} } 
    //      }
    //
    //      First bunch spans across 4 MPI processors of rank (0, 1, 
    //      2, and 3), in the context of the bunch simulator 
    //      communicator. The second bunch on rank (4, 5, 6, and 7). 
    //
    //      bunch_idx_map = { {0}, {-1} } on rank 0
    //      bunch_idx_map = { {0}, {-1} } on rank 1
    //      bunch_idx_map = { {0}, {-1} } on rank 2
    //      bunch_idx_map = { {0}, {-1} } on rank 3
    //      bunch_idx_map = { {-1}, {0} } on rank 4
    //      bunch_idx_map = { {-1}, {0} } on rank 5
    //      bunch_idx_map = { {-1}, {0} } on rank 6
    //      bunch_idx_map = { {-1}, {0} } on rank 7
    //
    //
    //   2. num_bunches_pri = 2, num_bunches_sec = 2
    //      4 bunches, 4 ranks (n_ranks = n_bunches)
    //
    //      bunch_ranks = { 
    //          { {0}, {1} }, 
    //          { {2}, {3} } 
    //      }
    //
    //      Each of the bunch takes over a single rank from 0 to 3. 
    //
    //      bunch_idx_map = { { 0, -1}, {-1, -1} } on rank 0
    //      bunch_idx_map = { {-1,  0}, {-1, -1} } on rank 1
    //      bunch_idx_map = { {-1, -1}, { 0, -1} } on rank 2
    //      bunch_idx_map = { {-1, -1}, {-1,  0} } on rank 3
    //
    //
    //   3. num_bunches_pri = 2, num_bunches_sec = 2
    //      4 bunches, 2 ranks (n_ranks < n_bunches)
    //
    //      bunch_ranks = { 
    //          { {0}, {0} }, 
    //          { {1], {1} }
    //      }
    //
    //      The first two bunches reside on rank 0, and the next two 
    //      on rank 1. 
    //
    //      bunch_idx_map = { { 0,  1}, {-1, -1} } on rank 0
    //      bunch_idx_map = { {-1, -1}, { 0,  1} } on rank 1
    //
    //   4. num_bunches_pri = 4, num_bunches_sec = 2
    //      6 bunches, 3 ranks (n_ranks < n_bunches)
    //
    //      bunch_ranks = { 
    //          { {0}, {0}, {1}, {1} }, 
    //          { {2}, {2} }
    //      }
    //
    //      bunch_idx_map = { { 0,  1, -1, -1}, {-1, -1} } on rank 0
    //      bunch_idx_map = { {-1, -1,  0,  1}, {-1, -1} } on rank 1
    //      bunch_idx_map = { {-1, -1, -1, -1}, { 0,  1} } on rank 2
    //
}


Bunch_simulator::Bunch_simulator(
        Bunch_train && pt, 
        Bunch_train && st,
        Commxx const & comm )
    : trains{std::move(pt), std::move(st)}
    , comm(comm)
    , diags_step_period()
    , diags_step_listed()
    , diags_element()
    , prop_actions()
    , prop_actions_step_end()
    , prop_actions_turn_end()
{
}

std::vector<int> 
Bunch_simulator::get_bunch_ranks(size_t train, size_t bunch) const
{
    int rank_per_bunch = std::ceil( 1.0 * comm.size() / 
        (trains[0].get_num_bunches() + trains[1].get_num_bunches()) );

    std::vector<int> ranks(rank_per_bunch);
    std::iota(ranks.begin(), ranks.end(), bunch*rank_per_bunch);

    return ranks;
}

void
Bunch_simulator::diag_action_step_and_turn(int turn_num, int step_num)
{
    for (auto const& dt : diags_step_period)
    {
        if (dt.trigger(turn_num, step_num))
        {
            get_bunch(dt.train, dt.bunch)
                .diag_update_and_write(dt.diag_name);
        }
    }

    for (auto const& dt : diags_step_listed)
    {
        if (dt.trigger(turn_num, step_num))
        {
            get_bunch(dt.train, dt.bunch)
                .diag_update_and_write(dt.diag_name);
        }
    }
}

void
Bunch_simulator::diag_action_element(Lattice_element const& element)
{
    for (auto const& dt : diags_element)
    {
        if (dt.trigger(element))
        {
            get_bunch(dt.train, dt.bunch)
                .diag_update_and_write(dt.diag_name);
        }
    }

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
Bunch_simulator::reg_prop_action_step_end(action_step_t fun)
{
    prop_actions_step_end.push_back(fun);
}

void
Bunch_simulator::reg_prop_action_turn_end(action_turn_t fun)
{
    prop_actions_turn_end.push_back(fun);
}

void
Bunch_simulator::reg_prop_action_step_end(action_data_step_t fun, void* data)
{
    using namespace std::placeholders;
    auto fun2 = std::bind(fun, _1, _2, _3, _4, data);
    prop_actions_step_end.push_back(fun2);
}

void
Bunch_simulator::reg_prop_action_turn_end(action_data_turn_t fun, void* data)
{
    using namespace std::placeholders;
    auto fun2 = std::bind(fun, _1, _2, _3, data);
    prop_actions_turn_end.push_back(fun2);
}

void
Bunch_simulator::prop_action_first(Lattice & lattice)
{
    if (prop_actions) 
        prop_actions->first(*this, lattice);
}

void
Bunch_simulator::prop_action_step_end(Lattice & lattice, int turn, int step)
{
    if (prop_actions) 
        prop_actions->step_end(*this, lattice, turn, step);

    for (auto const& action : prop_actions_step_end)
        action(*this, lattice, turn, step);
}

void
Bunch_simulator::prop_action_turn_end(Lattice & lattice, int turn)
{
    if (prop_actions) 
        prop_actions->turn_end(*this, lattice, turn);

    for (auto const& action : prop_actions_turn_end)
        action(*this, lattice, turn);
}

void 
Bunch_simulator::set_lattice_reference_particle(Reference_particle const & ref)
{
    for(auto & train : get_trains())
    {
        for(auto & bunch : train.get_bunches())
        {
            bunch.set_design_reference_particle(ref);
        }
    }
}
