#ifndef SPLIT_OPERATOR_STEPPER_CHOICE_H_
#define SPLIT_OPERATOR_STEPPER_CHOICE_H_
#include "synergia/simulation/stepper.h"

struct Kicks
{
    Kicks() :
                    num_steps(0)
    {
    }
    ;

    Kicks(Collective_operators const& collective_operators, int num_steps) :
                    collective_operators(collective_operators),
                    num_steps(num_steps)
    {
    }

    Collective_operators collective_operators;
    int num_steps;

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(Kicks);
typedef std::map<std::string, Kicks > List_choice_map;

/// Generate steps according with a list
class Split_operator_stepper_choice : public Stepper
{
private:
    List_choice_map list_choice_map;
    int num_steps_else;
    void
    make_stepper_else(Lattice_elements::iterator const& begin,
            Lattice_elements::iterator const& end,
            double const & length_between, double const & max_step_length);
    void
    construct_per_element_else();
    void
    construct_split_else();
public:
    Split_operator_stepper_choice(Lattice_sptr lattice_sptr, int map_order,
            List_choice_map const & list_choice_map, bool split_else = true);
    Split_operator_stepper_choice(int num_steps_else,
            Lattice_sptr lattice_sptr, int map_order,
            List_choice_map const & list_choice_map, bool split_else = true);
    /// Deprecated.
    Split_operator_stepper_choice(Lattice_simulator const& lattice_simulator,
            List_choice_map const & list_choice_map, bool split_else = true);
    /// Deprecated.
    Split_operator_stepper_choice(int num_steps_else,
            Lattice_simulator const& lattice_simulator,
            List_choice_map const & list_choice_map, bool split_else = true);

    /// Default constructor for serialization use only
    Split_operator_stepper_choice();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Split_operator_stepper_choice();
};
BOOST_CLASS_EXPORT_KEY(Split_operator_stepper_choice);
typedef boost::shared_ptr<Split_operator_stepper_choice > Split_operator_stepper_choice_sptr;

#endif /* SPLIT_OPERATOR_STEPPER_CHOICE_H_ */
