#ifndef SIMULATION_BUNCH_SIMULATOR_H
#define SIMULATION_BUNCH_SIMULATOR_H

#include "synergia/bunch/bunch_train.h"
#include "synergia/lattice/lattice.h"

class Operator;
class Independent_operation;

class Bunch_simulator
{

public:

    static const int FINAL_STEP = -1;

private:

    enum class diag_event_t
    {
        step_and_turn,
        element,
        syn_operator,
        syn_operation,
        particle_loss,
    };

    // bool trigger(int turn, int step)
    using trigger_step_t = std::function<bool(int, int)>;

    // bool trigger(void)
    using trigger_loss_t = std::function<bool(void)>;

    // bool trigger(Lattice_element const & ele)
    using trigger_ele_t  = std::function<bool(Lattice_element const &)>;

    // bool trigger(Operator const & opr)
    using trigger_opr_t  = std::function<bool(Operator const &)>;

    // bool trigger(Operation const & opn)
    using trigger_opn_t  = std::function<bool(Independent_operation const &)>;


    template<typename TriggerT>
    struct diag_tuple_t
    {
        int train;
        int bunch;
        std::string diag_name;
        TriggerT trigger;
    };

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

    // factory methods
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


    // diagnostics registration
    template<typename DiagT>
    void reg_diag(
            std::string const& name,
            DiagT const & diag, 
            trigger_step_t trig, 
            int train, 
            int bunch )
    { 
        get_bunch(train, bunch).add_diagnostics(name, diag);
        diags_step.emplace_back(train, bunch, name, trig);
    }

#if 0
    template<typename DiagT>
    void reg_diag_particle_loss(
            DiagT & diag, 
            trigger_loss_t,
            int train = 0, int bunch = 0 );

    template<typename DiagT>
    void reg_diag_element(
            DiagT & diag, 
            trigger_ele_t,
            int train = 0, int bunch = 0 );

    template<typename DiagT>
    void reg_diag_operator(
            DiagT & diag, 
            trigger_opr_t,
            int train = 0, int bunch = 0 );

    template<typename DiagT>
    void reg_diag_operation(
            DiagT & diag, 
            trigger_opn_t,
            int train = 0, int bunch = 0 );
#endif

    template<typename DiagT>
    void reg_diag_per_turn(
            std::string const& name,
            DiagT const & diag, 
            int train = 0, 
            int bunch = 0, 
            int period = 1 )
    { 
        auto trig = [period](int turn, int step) { 
            return step==Bunch_simulator::FINAL_STEP && turn%period==0; 
        };

        reg_diag( name, diag, trig, train, bunch );
    }

    void diag_action_step_and_turn(int turn_num, int step_num);
    void diag_action_particle_loss_update( );
    void diag_action_particle_loss_write( );
    void diag_action_element(Lattice_element const & element);
    void diag_action_operator(Operator const & opr);
    void diag_action_operation(Independent_operation const & opn);

    void prop_action_first(Lattice & lattice);
    void prop_action_step_end(Lattice & lattice, int turn, int step);
    void prop_action_turn_end(Lattice & lattice, int turn);

    void set_lattice_reference_particle(Reference_particle const & ref);

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

    std::vector<diag_tuple_t<trigger_step_t>> diags_step;
    std::vector<diag_tuple_t<trigger_loss_t>> diags_loss;
    std::vector<diag_tuple_t<trigger_ele_t>>  diags_ele;
    std::vector<diag_tuple_t<trigger_opr_t>>  diags_opr;
    std::vector<diag_tuple_t<trigger_opn_t>>  diags_opn;
};

#endif
