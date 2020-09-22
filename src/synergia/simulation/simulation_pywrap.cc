
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"

#include "synergia/simulation/stepper.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/simulation/checkpoint.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(simulation, m)
{
    m.def( "init", [](){ 
            //MPI_Init(NULL, NULL);
            Kokkos::initialize(); 
         });

    m.def( "finalize", [](){ 
            //MPI_Finalize(); 
            Kokkos::finalize(); 
         });

    m.def( "checkpoint_save",
            &syn::checkpoint_save,
            "propagator"_a,
            "bunch_simulator"_a )
        ;

    m.def( "checkpoint_load",
            &syn::checkpoint_load )
        ;
            
    // Propagator
    py::class_<Propagator>(m, "Propagator")
        .def(py::init<Lattice const&, Stepper const&>())
        .def("propagate", &Propagator::propagate)
        .def("print_steps", &Propagator::print_steps)
        ;

    // Lattice simulator -- only a namespace
    m.def_submodule("Lattice_simulator")
        .def( "tune_linear_lattice", 
                &Lattice_simulator::tune_linear_lattice,
                "Tune linear lattice.",
                "lattice"_a )

        .def( "tune_circular_lattice", 
                &Lattice_simulator::tune_circular_lattice,
                "Tune a circular lattice.",
                "lattice"_a, "tolerance"_a = 1.0e-13 )
        ;

    // Collective operator options (base class)
    py::class_<CO_options>(m, "CO_options");

    // Stepper base class
    py::class_<Stepper>(m, "Stepper")
        .def_readonly_static( "fixed_step_tolerance", 
                &Stepper::fixed_step_tolerance)
        ;

    // Split_operator_stepper
    py::class_<Split_operator_stepper, Stepper>(m, "Split_operator_stepper")
        .def( py::init<CO_options const&, int>(),
                "collective_options"_a, "num_steps"_a )
        .def( "append_collective_op",
                &Split_operator_stepper::append_collective_op,
                "Append a collective operator",
                "collective_options"_a )
        ;

    // Independent_stepper_elements
    py::class_<Independent_stepper_elements, Stepper>(m, "Independent_stepper_elements")
        .def( py::init<int>(), "steps_per_element"_a = 1 )
        ;


    // Bunch_simulator
    using action_step_t = std::function<void(Bunch_simulator&, Lattice&, int, int)>;
    using action_turn_t = std::function<void(Bunch_simulator&, Lattice&, int)>;

    py::class_<Bunch_simulator>(m, "Bunch_simulator")
        .def_static( "create_single_bunch_simulator",
                &Bunch_simulator::create_single_bunch_simulator,
                py::return_value_policy::move,
                "Create a Bunch_simulator with a single bunch.",
                "reference_particle"_a, 
                "num_particles"_a, 
                "num_real_particles"_a, 
                "commxx"_a = Commxx() )
 
        .def_static( "create_bunch_train_simulator",
                &Bunch_simulator::create_bunch_train_simulator,
                py::return_value_policy::move,
                "Create a Bunch_simulator with a single bunch train.",
                "reference_particle"_a, 
                "num_particles"_a, 
                "num_real_particles"_a, 
                "num_bunches"_a,
                "spacing"_a,
                "commxx"_a = Commxx() )
        
#if 0
        .def( "set_turns",
                &Bunch_simulator::set_turns,
                "Set the simulation start turn and total number of turns.",
                "first_turn"_a, "num_turns"_a )
#endif

        .def( "get_bunch",
                (Bunch const& (Bunch_simulator::*)(size_t, size_t) const)&Bunch_simulator::get_bunch,
                //py::overload_cast<size_t, size_t>(&Bunch_simulator::get_bunch, py::const_),
                py::return_value_policy::reference_internal,
                "Get the bunch reference from the bunch_simulator.",
                "train"_a = 0, "bunch"_a = 0 )

        .def( "reg_prop_action_step_end",
                (void (Bunch_simulator::*)(action_step_t))&Bunch_simulator::reg_prop_action_step_end,
                //py::overload_cast<action_step_t>(&Bunch_simulator::reg_prop_action_step_end),
                "Register the step end propagate action (callback)."
                "\nthe callback should have a signature of (Bunch_simulator, Lattice, turn_num, step_num)",
                "action"_a )

        .def( "reg_prop_action_turn_end",
                (void (Bunch_simulator::*)(action_turn_t))&Bunch_simulator::reg_prop_action_turn_end,
                //py::overload_cast<action_turn_t>(&Bunch_simulator::reg_prop_action_turn_end),
                "Register the turn end propagate action (callback)."
                "\nthe callback should have a signature of (Bunch_simulator, Lattice, turn_num)",
                "action"_a )

        .def( "reg_diag_per_turn", 
                &Bunch_simulator::reg_diag_per_turn<std::shared_ptr<Diagnostics>>,
                "Register a per turn diagnostics.",
                "diag"_a, 
                "train_idx"_a = 0, 
                "bunch_idx"_a = 0, 
                "period"_a = 1 )

        .def( "reg_diag_per_step", 
                &Bunch_simulator::reg_diag_per_step<std::shared_ptr<Diagnostics>>,
                "Register a per step diagnostics.",
                "diag"_a, 
                "train_idx"_a = 0, 
                "bunch_idx"_a = 0, 
                "period"_a = 1 )

        .def( "reg_diag_turn_listed", 
                &Bunch_simulator::reg_diag_turn_listed<std::shared_ptr<Diagnostics>>,
                "Register a per step diagnostics.",
                "diag"_a, 
                "train_idx"_a = 0, 
                "bunch_idx"_a = 0, 
                "turns"_a = std::vector<int>() )

        .def( "current_turn",
                &Bunch_simulator::current_turn,
                "Get the current simulation turn." )

        .def( "max_turns",
                &Bunch_simulator::max_turn,
                "Get the max simulation turns." )

        .def( "set_max_turns",
                &Bunch_simulator::set_max_turn,
                "Set the max simulation turns.",
                "max_turns"_a )

        .def( "dump",
                &Bunch_simulator::dump,
                "Dump." )

        ;


}


