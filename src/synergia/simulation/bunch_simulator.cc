#include <numeric>

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/independent_operation.h"

namespace impl {

    void divide_bunches(int size, 
            size_t num_bunches_pri, 
            size_t num_bunches_sec,
            std::vector<int> & p_ranks,
            std::vector<int> & s_ranks )
    {
        assert(size > 0);

        const size_t num_bunches = num_bunches_pri + num_bunches_sec;
        size_t st_start_rank = 0;

        if ( num_bunches == 0 )
        {
            // no bunch at all
            p_ranks.resize(0);
            s_ranks.resize(0);
            st_start_rank = 0;
        } 
        else if ( size == 1 )
        {
            // all bunches on a single rank
            p_ranks.resize(num_bunches_pri ? 1 : 0);
            s_ranks.resize(num_bunches_sec ? 1 : 0);
            st_start_rank = 0;
        }
        else if ( size < num_bunches )
        {
            // multiple ranks, each rank stores 1 train at max
            if (num_bunches % size != 0)
            {
                throw std::runtime_error(
                        "Bunch_simulator::create_bunch_train_simulator() "
                        "the number of bunches must be divisible by the number of ranks" );
            }

            int const bunch_per_rank = num_bunches / size;

            if (num_bunches_pri % bunch_per_rank != 0
                    || num_bunches_sec % bunch_per_rank != 0 )
            {
                throw std::runtime_error(
                        "Bunch_simulator::create_bunch_train_simulator() "
                        "the number of bunches in primary or secondary train must be "
                        "divisible by the number of bunches per rank" );
            }

            p_ranks.resize(num_bunches_pri/bunch_per_rank);
            s_ranks.resize(num_bunches_sec/bunch_per_rank);
            st_start_rank = p_ranks.size();
        }
        else
        {
            // now size >= num_bunches, one bunch (or a fraciton of a bunch) per rank
            if (size % num_bunches != 0)
            {
                throw std::runtime_error(
                        "Bunch_simulator::create_bunch_train_simulator() "
                        "the number of ranks must be divisible by the number of bunches" );
            }

            int const rank_per_bunch = size / num_bunches;

            p_ranks.resize(rank_per_bunch * num_bunches_pri);
            s_ranks.resize(rank_per_bunch * num_bunches_sec); 
            st_start_rank = p_ranks.size();
        }

        std::iota(p_ranks.begin(), p_ranks.end(), 0);
        std::iota(s_ranks.begin(), s_ranks.end(), st_start_rank);

        return;
    }

}

Bunch_simulator
Bunch_simulator::create_empty_bunch_simulator()
{
    return construct( Reference_particle(), 
                      Reference_particle(),
                      1.0, 1.0,  // num_particle, num_real_particle
                      0,   0,    // num bunches
                      1.0, 1.0,  // spacing
                      Commxx() );
}

Bunch_simulator
Bunch_simulator::create_single_bunch_simulator(
        Reference_particle const& ref,
        size_t num_particles,
        double num_real_particles,
        Commxx const& comm )
{
    return construct( ref, ref, 
                      num_particles, 
                      num_real_particles,
                      1,   0,   // num_bunches
                      1.0, 1.0, // spacing
                      comm );
}


Bunch_simulator 
Bunch_simulator::create_bunch_train_simulator(
        Reference_particle const& ref,
        size_t num_particles,
        double num_real_particles,
        size_t num_bunches,
        double spacing,
        Commxx const& comm )
{
    return construct( ref, ref, 
                      num_particles, 
                      num_real_particles,
                      num_bunches, 0, // num_bunches
                      spacing, 1.0,   // spacing
                      comm );
}

Bunch_simulator 
Bunch_simulator::create_two_trains_simulator(
        Reference_particle const& ref_pri,
        Reference_particle const& ref_sec,
        size_t num_particles,
        double num_real_particles,
        size_t num_bunches_pri,
        size_t num_bunches_sec,
        double spacing_pri,
        double spacing_sec,
        Commxx const& comm )
{
    return construct( ref_pri, ref_sec, 
                      num_particles, 
                      num_real_particles,
                      num_bunches_pri,
                      num_bunches_sec,
                      spacing_pri, 
                      spacing_sec,
                      comm );
}


Bunch_simulator
Bunch_simulator::construct(
        Reference_particle const& ref_pri,
        Reference_particle const& ref_sec,
        size_t num_part,
        double num_real_part,
        size_t num_bunches_pri,
        size_t num_bunches_sec,
        double spacing_pri,
        double spacing_sec,
        Commxx const& comm)
{
    auto comm_ptr = std::make_shared<Commxx>(comm);

    std::vector<int> p_ranks;
    std::vector<int> s_ranks;

    impl::divide_bunches( comm.size(),
            num_bunches_pri, 
            num_bunches_sec,
            p_ranks, 
            s_ranks );

    auto comm_pri = comm_ptr->group(p_ranks);
    auto comm_sec = comm_ptr->group(s_ranks);

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
                            comm_ptr );



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
        std::shared_ptr<Commxx> const& comm )
    : trains{std::move(pt), std::move(st)}
    , comm(std::move(comm))
    , diags_step_period()
    , diags_step_listed()
    , diags_element()
    , prop_actions()
    , prop_actions_step_end()
    , prop_actions_turn_end()
{
}

int 
Bunch_simulator::get_bunch_array_idx(int train, int bunch) const
{ 
    if (train > 1) return -1;
    if (bunch >= trains[train].get_num_bunches()) return -1;

    return trains[train].get_bunch_array_idx(bunch);
}

bool
Bunch_simulator::has_local_bunch(size_t train, size_t bunch) const
{
    return get_bunch_array_idx(train, bunch) != -1;
}

Bunch & 
Bunch_simulator::get_bunch(size_t train, size_t bunch)
{
    auto idx = get_bunch_array_idx(train, bunch);
    if (idx == -1) throw std::runtime_error("bunch not avaialble on the rank");
    return trains[train][idx];
}

Bunch const& 
Bunch_simulator::get_bunch(size_t train, size_t bunch) const
{
    auto idx = get_bunch_array_idx(train, bunch);
    if (idx == -1) throw std::runtime_error("bunch not avaialble on the rank");
    return trains[train][idx];
}

std::vector<int> 
Bunch_simulator::get_bunch_ranks(size_t train, size_t bunch) const
{
    int rank_per_bunch = std::ceil( 1.0 * comm->size() / 
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
            trains[dt.train][dt.bunch]
                .diag_update_and_write(dt.diag_name);
        }
    }

    for (auto const& dt : diags_step_listed)
    {
        if (dt.trigger(turn_num, step_num))
        {
            trains[dt.train][dt.bunch]
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
            trains[dt.train][dt.bunch]
                .diag_update_and_write(dt.diag_name);
        }
    }
}


void
Bunch_simulator::diag_action_operator(Operator const& opr)
{
}

void
Bunch_simulator::diag_action_operation(Independent_operation const& opn)
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
Bunch_simulator::set_lattice_reference_particle(Reference_particle const& ref)
{
    for(auto & train : get_trains())
    {
        for(auto & bunch : train.get_bunches())
        {
            bunch.set_design_reference_particle(ref);
        }
    }
}

namespace 
{
    std::vector<Bunch *> 
    get_bunch_ptrs(std::array<Bunch_train, 2> & trains)
    {
        std::vector<Bunch *> bunches;

        for(auto & t : trains)
            for(auto & b : t.get_bunches())
                bunches.push_back(&b);

        return bunches;
    }
}

void
Bunch_simulator::save_bunch_particles()
{
    Hdf5_file file("bunch_simulator.h5", Hdf5_file::truncate, *comm);
    auto bunches = get_bunch_ptrs(trains);

    for (int i=0; i<bunches.size(); ++i)
        bunches[i]->save_particles(file, i);
}

void
Bunch_simulator::load_bunch_particles()
{
    Hdf5_file file("bunch_simulator.h5", Hdf5_file::read_only, *comm);
    auto bunches = get_bunch_ptrs(trains);

    for (int i=0; i<bunches.size(); ++i)
        bunches[i]->load_particles(file, i);
}

