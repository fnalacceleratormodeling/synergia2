#ifndef SIMULATION_BUNCH_SIMULATOR_H
#define SIMULATION_BUNCH_SIMULATOR_H

#include "synergia/bunch/bunch_train.h"
#include "synergia/lattice/lattice.h"
#include "synergia/simulation/propagate_actions.h"

#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>  // std::pair

class Operator;
class Independent_operation;

class Diagnostics_loss;
class Bunch_simulator;

class Bunch_simulator
{

public:

    static const int FINAL_STEP = -1;

private:

#if 0
    // bool trigger(int turn, int step)
    using trigger_step_t = std::function<bool(int, int)>;

    // bool trigger(Lattice_element const & ele)
    using trigger_ele_t  = std::function<bool(Lattice_element const &)>;

    // bool trigger(Operator const & opr)
    using trigger_opr_t  = std::function<bool(Operator const &)>;

    // bool trigger(Operation const & opn)
    using trigger_opn_t  = std::function<bool(Independent_operation const &)>;
#endif

    struct trigger_step_period
    {
        int turn_period, step_period;

        bool operator()(int turn, int step) const
        {
            return (step == FINAL_STEP) && (turn % turn_period == 0); 
        }

        template<class AR>
        void serialize(AR & ar)
        { ar(turn_period, step_period); }
    };

    struct trigger_step_listed
    {
        std::vector<std::pair<int, int>> numbers;

        bool operator()(int turn, int step) const
        {
            return false;
        }

        template<class AR>
        void serialize(AR & ar)
        { ar(numbers); }
    };

    struct trigger_element
    {
        std::vector<std::string> names;

        bool operator()(Lattice_element const& ele) const
        {
            return false;
        }

        template<class AR>
        void serialize(AR & ar)
        { ar(names); }
    };

    template<typename TriggerT>
    struct diag_tuple_t
    {
        int train;
        int bunch;
        std::string diag_name;
        TriggerT trigger;

        template<class AR>
        void serialize(AR & ar)
        { ar(train, bunch, diag_name, trigger); }
    };

    using dt_step_period = diag_tuple_t<trigger_step_period>;
    using dt_step_listed = diag_tuple_t<trigger_step_listed>;
    using dt_element     = diag_tuple_t<trigger_element>;

    // void action(Bunch_simulator&, Lattice&, int turn, int step, void* data)
    using action_step_t = std::function<void(Bunch_simulator&, Lattice&, int, int)>;
    using action_data_step_t = std::function<void(Bunch_simulator&, Lattice&, int, int, void*)>;

    // void action(Bunch_simulator&, Lattice&, int turn, void* data)
    using action_turn_t = std::function<void(Bunch_simulator&, Lattice&, int)>;
    using action_data_turn_t = std::function<void(Bunch_simulator&, Lattice&, int, void*)>;

private:

    // constructor
    Bunch_simulator( Bunch_train && pt, 
                     Bunch_train && st, 
                     std::vector<std::vector<int>> const & bunch_ranks,
                     Commxx const & comm );


    static Bunch_simulator 
        construct( Reference_particle const & ref_pri,
                   Reference_particle const & ref_sec,
                   size_t num_part,
                   double num_real_part,
                   size_t num_bunches_pri,
                   size_t num_bunches_sec,
                   double spacing_pri,
                   double spacing_sec,
                   Commxx const & comm);



public:

    Bunch_simulator(Bunch_simulator const&) = delete;
    Bunch_simulator(Bunch_simulator &&) = default;

    // factory methods
    static Bunch_simulator 
        create_empty_bunch_simulator();

    static Bunch_simulator 
        create_single_bunch_simulator(
            Reference_particle const & ref,
            size_t num_particles,
            double num_real_particles,
            Commxx const & comm = Commxx()
            );

    static Bunch_simulator 
        create_bunch_train_simulator(
            Reference_particle const & ref,
            size_t num_particles,
            double num_real_particles,
            size_t num_bunches,
            double spacing,
            Commxx const & comm = Commxx()
            );

    static Bunch_simulator 
        create_two_trains_simulator(
            Reference_particle const & primary_ref,
            Reference_particle const & secondary_ref,
            size_t num_particles,
            double num_real_particles,
            size_t num_bunches_primary,
            size_t num_bunches_secondary,
            Commxx const & comm = Commxx()
            );


    // diag per turn
    template<class Diag>
    void reg_diag_per_turn(
            Diag const& diag, 
            std::string const& name, 
            std::string const& filename,
            int train = 0, int bunch = 0, int period = 1 )
    { 
        dt_step_period dt{train, bunch, name, trigger_step_period{period, -1}};
        diags_step_period.push_back(dt);
        get_bunch(train, bunch).add_diagnostics(diag, name, filename);
    }

#if 0
    // diag loss
    void reg_diag_loss_aperture(
            Diagnostics_loss & diag, int train = 0, int bunch = 0 )
    { get_bunch(train, bunch).set_diag_loss_aperture(diag); }

    void reg_diag_loss_zcut(
            Diagnostics_loss & diag, int train = 0, int bunch = 0 )
    { get_bunch(train, bunch).set_diag_loss_zcut(diag); }

    // diag per element
    void reg_diag_element(
            std::string const& name, Diagnostics & diag, 
            trigger_ele_t, int train = 0, int bunch = 0 )
    { }

    // diag per operator
    void reg_diag_operator(
            std::string const& name, Diagnostics & diag, 
            trigger_opr_t, int train = 0, int bunch = 0 )
    { }

    // diag per operation
    void reg_diag_operation(
            std::string const& name, Diagnostics & diag, 
            trigger_opn_t, int train = 0, int bunch = 0 )
    { }
#endif

    // register propagation actions
    template<class PA>
    void reg_prop_actions(PA const& pa)
    { prop_actions = std::make_shared<PA>(pa); }

    // register individual propagation actions
    void reg_prop_action_step_end(action_step_t fun);
    void reg_prop_action_turn_end(action_turn_t fun);

    void reg_prop_action_step_end(action_data_step_t fun, void* data = nullptr);
    void reg_prop_action_turn_end(action_data_turn_t fun, void* data = nullptr);

    // diag actions
    void diag_action_step_and_turn(int turn_num, int step_num);
    void diag_action_element(Lattice_element const & element);
    void diag_action_operator(Operator const & opr);
    void diag_action_operation(Independent_operation const & opn);

    // propagate actions
    void prop_action_first(Lattice & lattice);
    void prop_action_step_end(Lattice & lattice, int turn, int step);
    void prop_action_turn_end(Lattice & lattice, int turn);

    // accessors
    std::array<Bunch_train, 2> & get_trains()
    { return trains; }

    std::array<Bunch_train, 2> const & get_trains() const
    { return trains; }

    Bunch_train & operator[](size_t idx)
    { return trains[idx]; }

    Bunch_train const & operator[](size_t idx) const
    { return trains[idx]; }

    Bunch & get_bunch(size_t train = 0, size_t bunch = 0)
    { return trains[train][bunch]; }

    Bunch const & get_bunch(size_t train = 0, size_t bunch = 0) const
    { return trains[train][bunch]; }

    std::vector<int> const & get_bunch_ranks(size_t train, size_t bunch) const
    { return bunch_ranks[train==0 ? bunch : pt_bunches + bunch]; }

    void set_turns(int first, int turns)
    { first_turn = first; num_turns = turns; }

    void set_lattice_reference_particle(Reference_particle const & ref);

public:

    int num_turns = 0;
    int first_turn = 0;
    int max_turns = 0;

private:

    std::array<Bunch_train, 2> trains;
    std::vector<std::vector<int>> bunch_ranks;
    int pt_bunches;
    int st_bunches;
    Commxx comm;

    std::vector<dt_step_period> diags_step_period;
    std::vector<dt_step_listed> diags_step_listed;
    std::vector<dt_element>     diags_element;

    // it would survive the checkpoint load
    std::shared_ptr<Propagate_actions> prop_actions;

    // non-persistent propagate actions -- 
    // reg again after checkpoint load
    std::vector<action_step_t> prop_actions_step_end;
    std::vector<action_turn_t> prop_actions_turn_end;

    friend class cereal::access;

    template<class AR>
    void serialize(AR & ar)
    {
        ar(CEREAL_NVP(trains));
        ar(CEREAL_NVP(bunch_ranks));
        ar(CEREAL_NVP(pt_bunches));
        ar(CEREAL_NVP(st_bunches));
        ar(CEREAL_NVP(comm));

        ar(CEREAL_NVP(diags_step_period));
        ar(CEREAL_NVP(diags_step_listed));
        ar(CEREAL_NVP(diags_element));

        ar(CEREAL_NVP(prop_actions));
    }
};

template<>
inline
void Bunch_simulator::reg_prop_actions(std::shared_ptr<Propagate_actions> const& pa)
{ prop_actions = pa; }

#endif
