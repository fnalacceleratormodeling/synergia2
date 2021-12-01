#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/lattice/lattice.h"

//#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/utils/cereal.h"

#include <cereal/types/memory.hpp>

class Bunch_simulator;

class Propagator
{
public:

    static const int PRE_TURN   = -1;
    static const int FINAL_STEP = -1;

    struct Slice_iterator
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = int;
        using value_type = Lattice_element_slice;
        using pointer    = Lattice_element_slice*;
        using reference  = Lattice_element_slice&;

        Slice_iterator(Propagator& p) : prop(p) { }

        reference operator*() const { return *slice_it; }
        pointer operator->() { return &(*slice_it); }

        void first_slice()
        {
            step_it = prop.steps.begin();
            step_end = prop.steps.end();

            if (step_it == prop.steps.end()) return;

            opr_it = step_it->operators.begin();
            opr_end = step_it->operators.end();

            if(is_valid_opr())
            {
                auto o = dynamic_cast<Independent_operator*>(opr_it->get());
                slice_it = o->slices.begin();
                slice_end = o->slices.end();
                return;
            }

            advance_to_next_slice();
        }

        void end()
        {
            step_it = prop.steps.end();
        }

        Slice_iterator& operator++()
        {
            ++slice_it;

            auto o = dynamic_cast<Independent_operator*>(opr_it->get());
            if(slice_it == o->slices.end())
                advance_to_next_slice();

            return *this;
        }

        void advance_to_next_slice()
        {
            while(step_it != prop.steps.end())
            {
                if (opr_it == step_it->operators.end())
                {
                    ++step_it;
                    if (step_it == prop.steps.end()) return;

                    opr_it = step_it->operators.begin();
                    opr_end = step_it->operators.end();
                }
                else
                {
                    ++opr_it;
                }

                if (is_valid_opr())
                {
                    auto o = dynamic_cast<Independent_operator*>(opr_it->get());
                    slice_it = o->slices.begin();
                    slice_end = o->slices.end();
                    return;
                }
            }

            // the step is now at the end if reached here
        }

        bool is_valid_opr()
        {
            if (opr_it == step_it->operators.end()) return false;

            auto o = dynamic_cast<Independent_operator*>(opr_it->get());
            if (o && o->slices.size()) return true;

            return false;
        }

        // iterator++
        Slice_iterator operator++(int)
        {
            Slice_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool 
        operator==(const Slice_iterator& a, const Slice_iterator& b)
        { 
            if (&a.prop != &b.prop) 
                return false;

            if (a.step_it != b.step_it) 
                return false;

            if (a.step_it != a.step_end &&
                b.step_it != b.step_end &&
                a.opr_it  != b.opr_it)
                return false;

            if (a.step_it != a.step_end &&
                b.step_it != b.step_end &&
                a.opr_it  != a.opr_end &&
                b.opr_it  != b.opr_end &&
                a.slice_it!= b.slice_it)
                return false;

            return true;
        }

        friend bool 
        operator!=(const Slice_iterator& a, const Slice_iterator& b)
        { return !(a==b); }

    private:

        Propagator& prop;

        std::vector<Step>::iterator step_it;
        std::vector<Step>::iterator step_end;

        std::vector<std::shared_ptr<Operator>>::iterator opr_it;
        std::vector<std::shared_ptr<Operator>>::iterator opr_end;

        std::vector<Lattice_element_slice>::iterator slice_it;
        std::vector<Lattice_element_slice>::iterator slice_end;
    };


    struct Lattice_element_slices
    { 
        Lattice_element_slices(Propagator& p) 
            : prop(p) { }

        Slice_iterator begin()
        {
            Slice_iterator it(prop);
            it.first_slice();
            return it;
        }

        Slice_iterator end()
        { 
            Slice_iterator it(prop);
            it.end();
            return it;
        }

        Propagator& prop;
    };

private:

    Lattice lattice;
    std::vector<Step> steps;
    Lattice_element_slices slices;
    
    std::unique_ptr<Stepper> stepper_ptr;

private:

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

public:

    // given lattice and stepper
    Propagator(
            Lattice const& lattice, 
            Stepper const& stepper = Independent_stepper_elements(1) )
        : lattice(lattice)
        , steps()
        , slices(*this)
        , stepper_ptr(stepper.clone())
    {
        this->lattice.update();
        steps = stepper_ptr->apply(this->lattice);
    }

    // max_turns: number of turns in this propagate. -1 run to the end
    void propagate(Bunch_simulator & simulator, Logger & logger, int max_turns = -1);

    Lattice_element_slices&
    get_lattice_element_slices()
    { return slices; }

    std::list<Lattice_element> const&
    get_lattice_elements()
    { return lattice.get_elements(); }

    Lattice&       get_lattice()       { return lattice; }
    Lattice const& get_lattice() const { return lattice; }

    // print
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
    Propagator() : lattice(), steps(), slices(*this), stepper_ptr() { }

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
