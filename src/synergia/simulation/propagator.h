#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/lattice/lattice.h"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/utils/cereal.h"

#include <cereal/types/memory.hpp>

class Propagator
{
public:

    static const int FINAL_STEP = -1;

private:

    Lattice lattice;
    std::vector<Step> steps;
    
    std::unique_ptr<Stepper> stepper_ptr;

private:

	void construct();

	void do_before_start(
            Bunch_simulator & simulator, 
            Logger & logger);

	void do_start_repetition(
            Bunch_simulator & simulator);

	void do_turn_end(
            Bunch_simulator & simulator, 
            int turn_count, 
            Logger & logger);

	void do_step(
            Bunch_simulator & simulator,
            Step & step, 
            int step_count, 
            int turn_count, 
            Logger & logger);

	bool check_out_of_particles(
            Bunch_simulator const & simulator, 
            Logger & logger);

#if 0
	void
	checkpoint(State & state, Logger & logger, double & t);
#endif

public:

    // given lattice and stepper
    Propagator(
            Lattice const& lattice, 
            Stepper const& stepper = Independent_stepper_elements(1) )
        : lattice(lattice)
        , steps()
        , stepper_ptr(stepper.clone())
    {
        this->lattice.update();
        steps = stepper_ptr->apply(this->lattice);
    }

    void propagate(Bunch_simulator & simulator, Logger & logger);

    void print_steps(Logger & logger) const
    { for(auto const & s : steps) s.print(logger); }

    // dump to a string
    std::string dump() const
    {
        std::stringstream ss;

        {
            cereal::JSONOutputArchive ar(ss);
            ar(*this);
        }

        return ss.str();
    }

    // static method to load from a serialize string to void
    // the public default constructor
    static Propagator load_from_string(std::string const& str)
    {
        std::stringstream ss(str);
        cereal::JSONInputArchive ar(ss);

        Propagator p;
        ar(p);

        return p;
    }

private:

    // default ctor for serialization only
    Propagator() : lattice(), steps(), stepper_ptr() { }

    friend class cereal::access;

    template<class AR>
    void save(AR & ar) const
    {
        ar(CEREAL_NVP(lattice));
        ar(CEREAL_NVP(stepper_ptr));
    }

    template<class AR>
    void load(AR & ar)
    {
        ar(CEREAL_NVP(lattice));
        ar(CEREAL_NVP(stepper_ptr));

        lattice.update();
        steps = stepper_ptr->apply(lattice);
    }

#if 0
    Stepper_sptr
    get_stepper_sptr();

    void
    set_checkpoint_period(int period);

    int
    get_checkpoint_period() const;

    void
    set_checkpoint_dir(std::string const& directory_name);

    std::string const&
    get_checkpoint_dir() const;

    void
    set_checkpoint_with_xml(bool with_xml);

    bool
    get_checkpoint_with_xml() const;

    void
    set_final_checkpoint(bool final_checkpoint);

    bool
    get_final_checkpoint() const;

    void
    set_concurrent_io(int max);

    int
    get_concurrent_io() const;

    void
    set_num_threads(int nt);

    int
    get_num_threads() const;

    void
    propagate(State & state);

    /// jfa note: the lifetime of the pointers in state must
    ///           be managed manually
    State
    get_resume_state(std::string const& checkpoint_dir);

    void
    resume( std::string const& checkpoint_dir, 
            bool new_num_turns, int num_turns, 
            bool new_max_turns, int max_turns, 
            bool new_verbosity, int verbosity );

    void
    propagate(
            Bunch_simulator & bunch_simulator, 
            int num_turns,
            int max_turns = 0, 
            int verbosity = 1 );

    void
    propagate(
            Bunch_simulator & bunch_simulator,
            Propagate_actions & general_actions, 
            int num_turns,
            int max_turns = 0, 
            int verbosity = 1 );
#endif

#if 0
    template<class Archive>
    void
    serialize(Archive & ar, const unsigned int version);
#endif
};

#endif /* PROPAGATOR_H_ */
