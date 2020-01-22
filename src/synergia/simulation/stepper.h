#ifndef STEPPER_H_
#define STEPPER_H_

#include <cereal/types/polymorphic.hpp>
#include <cereal/access.hpp>

#include "synergia/simulation/step.h"

class Stepper
{
public:

    static const double fixed_step_tolerance;

public:

    Stepper();

    virtual ~Stepper() = default;
    virtual std::unique_ptr<Stepper> clone() const = 0;

    std::vector<Step> 
    apply(Lattice const & lattice) const
    { 
        auto steps = apply_impl(lattice); 
        create_operations(lattice, steps);
        return steps;
    }

private:

    void create_operations(
            Lattice const & lattice,
            std::vector<Step> & steps) const
    { 
        for(auto & step : steps) 
            step.create_operations(lattice); 
    }

    virtual std::vector<Step> 
    apply_impl(Lattice const & lattice) const = 0;

    friend class cereal::access;

    template<class Archive>
    void serialize(Archive & ar)
    { }

#if 0
protected:

    Independent_operator_sptr get_fixed_step(
            std::string const& name, 
            Lattice_elements::iterator & lattice_it, 
            double & left, 
            Lattice_elements::iterator const & lattice_end, 
            const double step_length, 
            double & offset_fudge, 
            bool end_on_force_diagnostics);

    Lattice_element_slices extract_slices(Steps const& steps);
#endif

#if 0
    /// Deprecated
    Stepper(Lattice_simulator const& lattice_simulator);

    Lattice_simulator & get_lattice_simulator();

    Steps & get_steps();


    void force_update_operations_no_collective();

    virtual void print() const;

    virtual void cs_step_lattice_functions();
    virtual void print_cs_step_betas();

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif

};

// include the archive types before registering the derived class
#include <cereal/archives/json.hpp>

// derived class type will be registered in the header of each
// concrete stepper types, e.g.,
// CEREAL_REGISTER_TYPE(Independent_stepper)


#endif /* STEPPER_H_ */
