#ifndef STEPPER_H_
#define STEPPER_H_

#include "synergia/utils/cereal.h"
#include "synergia/simulation/step.h"

class Stepper
{
public:

    static const std::string force_diagnostics_attribute;
    static const double fixed_step_tolerance;

private:


protected:

#if 0
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

public:

    Stepper();
    virtual ~Stepper() = default;

    virtual std::vector<Step> apply(Lattice const & lattice) const = 0;

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

#endif /* STEPPER_H_ */
