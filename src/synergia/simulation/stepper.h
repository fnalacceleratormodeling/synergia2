#ifndef STEPPER_H_
#define STEPPER_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/utils/serialization.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"

class Stepper
{
private:
    Lattice_simulator lattice_simulator;
    Steps steps;

protected:    
    Independent_operator_sptr
    get_fixed_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double step_length, double & offset_fudge);

public:
    Stepper(Lattice_simulator const& lattice_simulator);
    /// Default constructor for serialization use only
    Stepper();
    Lattice_simulator &
    get_lattice_simulator();
    Steps &
    get_steps();
    virtual void
    print() const;

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(steps);
            ar & BOOST_SERIALIZATION_NVP(lattice_simulator);
        }

    virtual
    ~Stepper();
};

typedef boost::shared_ptr<Stepper > Stepper_sptr;

/// The Independent_stepper class generates evenly-spaced Independent_operator
/// steps through a Lattice. No collective effects are included.
class Independent_stepper : public Stepper
{
public:
    /// Construct an Independent_stepper
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param num_steps the number of steps to take in the Lattice
    Independent_stepper(Lattice_simulator const& lattice_simulator,
            int num_steps);

    /// Default constructor for serialization use only
    Independent_stepper();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
        }

    virtual
    ~Independent_stepper();

};

/// The Independent_stepper_elements class generates a constant number of
/// Independent_operator steps per thick element. Thin elements are assigned
/// a single step each. No collective effects are included.
class Independent_stepper_elements : public Stepper
{
public:
    /// Construct an Independent_stepper
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param steps_per_element the number of steps per thick element
    Independent_stepper_elements(Lattice_simulator const& lattice_simulator,
            int steps_per_element);

    /// Default constructor for serialization use only
    Independent_stepper_elements();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
        }

    virtual
    ~Independent_stepper_elements();
};

typedef boost::shared_ptr<Independent_stepper_elements >
        Independent_stepper_elements_sptr;

/// The Split_operator_stepper class generates evenly-spaced split-operator
/// steps through a Lattice. One or more collective effects are included per
/// step.
class Split_operator_stepper : public Stepper
{
    void
    construct(Collective_operators const & collective_operators, int num_steps);
public:
    /// Construct a Split_operator_stepper with a single Collective_operator
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operator the Collective_operator to apply in each step
    /// @param num_steps the number of steps to take in the Lattice
    Split_operator_stepper(Lattice_simulator const& lattice_simulator,
            Collective_operator_sptr collective_operator, int num_steps);

    /// Construct a Split_operator_stepper with multiple Collective_operators
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operators the set of Collective_operators to apply in each step
    /// @param num_steps the number of steps to take in the Lattice
    Split_operator_stepper(Lattice_simulator const& lattice_simulator,
            Collective_operators const & collective_operators, int num_steps);

    /// Default constructor for serialization use only
    Split_operator_stepper();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
        }

    virtual
    ~Split_operator_stepper();
};

typedef boost::shared_ptr<Split_operator_stepper > Split_operator_stepper_sptr;

/// The Split_operator_stepper_elements class generates a constant number of
/// split-operator steps per thick element. Thin elements are assigned
/// a single step each. One or more collective effects are included per
/// step.
class Split_operator_stepper_elements : public Stepper
{
private:
    void
    construct(Collective_operators const & collective_operators,
            int steps_per_element);
public:
    /// Construct a Split_operator_stepper_elements with a single Collective_operator
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operator the Collective_operator to apply in each step
    /// @param steps_per_element the number of steps per thick element
            Split_operator_stepper_elements(
                    Lattice_simulator const& lattice_simulator,
                    Collective_operator_sptr collective_operator,
                    int steps_per_element);

    /// Construct a Split_operator_stepper_elements with multiple Collective_operators
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operators the set of Collective_operators to apply
    ///        in each step
    /// @param steps_per_element the number of steps per thick element
    Split_operator_stepper_elements(Lattice_simulator const& lattice_simulator,
            Collective_operators const & collective_operators,
            int steps_per_element);

    /// Default constructor for serialization use only
    Split_operator_stepper_elements();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
        }

    virtual
    ~Split_operator_stepper_elements();
};

/// Generate steps through lattice based on envelope shape.
/// Includes collective effects.
//class Split_operator_stepper_smart : public Stepper
//{
//
//};

struct Kicks
{
    
    Kicks():
    num_steps(0)    
    {
    }; 
    
    Kicks(Collective_operators const& collective_operators, int num_steps):
    collective_operators(collective_operators), num_steps(num_steps)
    {
    } 
    
    Collective_operators collective_operators; 
    int num_steps;
};

typedef  std::map<std::string, Kicks >  List_choice_map;

/// Generate steps according with a list 
class Split_operator_stepper_choice: public Stepper
{
private:
   List_choice_map list_choice_map; 
   int num_steps_else;   
   void                    
   make_stepper_else(Lattice_elements::iterator const& begin, Lattice_elements::iterator const& end, double const & length_between, 
                double const & max_step_length);
   void
   construct_per_element_else();
   void
   construct_split_else();
public:  
  
  Split_operator_stepper_choice (Lattice_simulator const& lattice_simulator, List_choice_map const & list_choice_map, bool split_else=true);
  Split_operator_stepper_choice (int num_steps_else, Lattice_simulator const& lattice_simulator, List_choice_map const & list_choice_map,  bool split_else=true);
 
 
 virtual
    ~Split_operator_stepper_choice(); 
};

#endif /* STEPPER_H_ */
