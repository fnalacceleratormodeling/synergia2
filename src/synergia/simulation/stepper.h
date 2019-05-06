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
public:
    static const std::string force_diagnostics_attribute;
    static const double fixed_step_tolerance;

private:
    Lattice_simulator lattice_simulator;
    Steps steps;
    bool have_step_betas;

protected:
    Independent_operator_sptr
    get_fixed_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double step_length, double & offset_fudge,
        bool end_on_force_diagnostics);
    Lattice_element_slices
    extract_slices(Steps const& steps);

public:
    Stepper(Lattice_sptr lattice_sptr, int map_order);
    /// Deprecated
    Stepper(Lattice_simulator const& lattice_simulator);
    /// Default constructor for serialization use only
    Stepper();
    Lattice_simulator &
    get_lattice_simulator();
    Steps &
    get_steps();
    void
    force_update_operations_no_collective();

    virtual void
    print() const;

    virtual void
    cs_step_lattice_functions();

    virtual void 
    print_cs_step_betas();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Stepper();
};
typedef boost::shared_ptr<Stepper > Stepper_sptr; // syndoc:include

#endif /* STEPPER_H_ */
