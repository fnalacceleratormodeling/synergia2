#include "operator.h"
#include "lattice_simulator.h"
#include "stepper.h"
#include "propagator.h"
#include "propagate_actions.h"
#include "standard_diagnostics_actions.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

struct Propagate_actions_callback : Propagate_actions
{
    Propagate_actions_callback(PyObject *p) :
        Propagate_actions(), self(p)
    {
    }
    Propagate_actions_callback(PyObject *p, const Propagate_actions& x) :
        Propagate_actions(x), self(p)
    {
    }
    void
    first_action(Stepper & stepper, Bunch & bunch)
    {
        call_method<void > (self, "first_action", stepper, bunch);
    }
    static void
    default_first_action(Propagate_actions& self_, Stepper & stepper,
            Bunch & bunch)
    {
        self_.Propagate_actions::first_action(stepper, bunch);
    }
    void
    turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num)
    {
        call_method<void > (self, "turn_end_action", stepper, bunch, turn_num);
    }
    static void
    default_turn_end_action(Propagate_actions& self_, Stepper & stepper,
            Bunch & bunch, int turn_num)
    {
        self_.Propagate_actions::turn_end_action(stepper, bunch, turn_num);
    }
    void
    step_end_action(Stepper & stepper, Step & step, Bunch & bunch,
            int turn_num, int step_num)
    {
        call_method<void > (self, "step_end_action", stepper, step, bunch,
                turn_num, step_num);
    }
    static void
    default_step_end_action(Propagate_actions& self_, Stepper & stepper,
            Step & step, Bunch & bunch, int turn_num, int step_num)
    {
        self_.Propagate_actions::step_end_action(stepper, step, bunch,
                turn_num, step_num);
    }
private:
    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads23,
        Propagator::propagate, 2, 3);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads34,
        Propagator::propagate, 3, 4);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads45,
        Propagator::propagate, 4, 5);



BOOST_PYTHON_MODULE(simulation)
{
    class_<Operator, Operator_sptr, boost::noncopyable >("Operator", no_init)
        .def("get_name", &Operator::get_name,
                return_value_policy<copy_const_reference >())
        .def("get_type", &Operator::get_type,
                return_value_policy<copy_const_reference >())
        .def("apply", &Operator::apply)
        .def("print_", &Operator::print)
        ;

    class_<Collective_operator, Collective_operator_sptr,
        boost::noncopyable, bases<Operator > >
        ("Collective_operator", no_init)
//        .def("get_name", &Collective_operator::get_name)
//        .def("get_type", &Collective_operator::get_type)
//        .def("apply", &Collective_operator::apply)
//        .def("print_", &Collective_operator::print)
        ;

    class_<Dummy_collective_operator, Dummy_collective_operator_sptr,
        bases<Collective_operator> >("Dummy_collective_operator",
                init<std::string const& >())
//        .def("get_name", &Collective_operator::get_name)
//        .def("get_type", &Collective_operator::get_type)
//        .def("apply", &Collective_operator::apply)
//        .def("print_", &Collective_operator::print)
        ;

    class_<Independent_operator, Independent_operator_sptr,
        bases<Operator > >("Independent_operator", init<std::string const&,
                Operation_extractor_map_sptr>())
//        .def("get_name", &Collective_operator::get_name)
//        .def("get_type", &Collective_operator::get_type)
//        .def("apply", &Collective_operator::apply)
//        .def("print_", &Collective_operator::print)
        .def("append_slice", &Independent_operator::append_slice)
        .def("get_slices", &Independent_operator::get_slices,
                return_value_policy<copy_const_reference >())
        ;

    class_<Lattice_simulator >("Lattice_simulator",
            init<Lattice_sptr, int >())
        .def("set_slices",
                &Lattice_simulator::set_slices)
        .def("get_map_order", &Lattice_simulator::get_map_order)
        .def("get_operation_extractor_map",
                &Lattice_simulator::get_operation_extractor_map_sptr)
        .def("get_lattice", &Lattice_simulator::get_lattice_sptr)
        .def("get_chef_lattice", &Lattice_simulator::get_chef_lattice_sptr)
        .def("get_bucket_length", &Lattice_simulator::get_bucket_length)
        .def("get_number_buckets",&Lattice_simulator::get_number_buckets)
        .def("update", &Lattice_simulator::update)
        ;



    void (Step::*apply1)(Bunch &) = &Step::apply;
    class_<Step, Step_sptr >("Step", init<double >())
//            .def("append", remember how to overload methods...)
            .def("apply",apply1)
            .def("get_operators",&Step::get_operators,
                    return_value_policy<copy_const_reference >())
            ;
    to_python_converter<Steps,
             container_conversions::to_tuple<Steps > >();

    class_<Stepper >("Stepper",init<Lattice_simulator const& >())
        .def("get_lattice_simulator", &Stepper::get_lattice_simulator,
                return_internal_reference< >())
        .def("get_steps", &Stepper::get_steps,
                return_value_policy<copy_non_const_reference >())
        .def("print_", &Stepper::print)
        ;

    class_<Split_operator_stepper, bases<Stepper > >("Split_operator_stepper",
            init<Lattice_simulator const&, Collective_operator_sptr, int >());

    class_<Split_operator_stepper_elements, bases<Stepper > >("Split_operator_stepper_elements",
            init<Lattice_simulator const&, Collective_operator_sptr, int >());

    class_<Independent_stepper, bases<Stepper > >("Independent_stepper",
            init<Lattice_simulator const&, int >());

    class_<Independent_stepper_elements, bases<Stepper > >("Independent_stepper_elements",
            init<Lattice_simulator const&, int >());

    class_<Propagate_actions, Propagate_actions_callback >("Propagate_actions",
            init< >())
            .def("first_action",
                    &Propagate_actions_callback::default_first_action)
            .def("turn_end_action",
                    &Propagate_actions_callback::default_turn_end_action)
            .def("step_end_action",
                    &Propagate_actions_callback::default_step_end_action)
            ;

    class_<Standard_diagnostics_actions, bases<Propagate_actions > >(
            "Standard_diagnostics_actions", init< >())
            .def("add_per_turn", &Standard_diagnostics_actions::add_per_turn)
            .def("add_per_step", &Standard_diagnostics_actions::add_per_step)
            .def("update_and_write_all", &Standard_diagnostics_actions::update_and_write_all)
            .def("first_action", &Standard_diagnostics_actions::first_action)
            .def("turn_end_action", &Standard_diagnostics_actions::turn_end_action)
            .def("step_end_action", &Standard_diagnostics_actions::step_end_action)
            ;

    void (Propagator::*propagate1)(Bunch_with_diagnostics &, int, bool) 
                                = &Propagator::propagate;

    void (Propagator::*propagate2)(Bunch_with_diagnostics_train &, int, bool) 
                                = &Propagator::propagate;                                

    void (Propagator::*propagate3)(Bunch &, int, Diagnostics &,
            Diagnostics &, bool) = &Propagator::propagate;

    void (Propagator::*propagate4)(Bunch &, int, Multi_diagnostics &,
            Multi_diagnostics &, bool) = &Propagator::propagate;
            
    void (Propagator::*propagate5)(Bunch &, int, Propagate_actions &,
            int) = &Propagator::propagate; 
               
    void (Propagator::*propagate6)(Bunch &, int, Propagate_actions &,
            Propagate_actions &, int) = &Propagator::propagate;               
                                
    class_<Propagator >("Propagator",init<Stepper_sptr >())
            .def("propagate", propagate1,
                 propagate_member_overloads23())
            .def("propagate", propagate2,
                    propagate_member_overloads23())
            .def("propagate", propagate3,
                    propagate_member_overloads45())
            .def("propagate", propagate4,
                    propagate_member_overloads45())
             .def("propagate", propagate5,
                    propagate_member_overloads34())
            .def("propagate", propagate6,
                    propagate_member_overloads45())          
            ;

}
