#include "operator.h"
#include "lattice_simulator.h"
#include "stepper.h"
#include "propagator.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads,
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
        .def("construct_sliced_chef_beamline",
                &Lattice_simulator::construct_sliced_chef_beamline)
        .def("get_map_order", &Lattice_simulator::get_map_order)
        .def("get_operation_extractor_map",
                &Lattice_simulator::get_operation_extractor_map_sptr)
        .def("get_lattice", &Lattice_simulator::get_lattice_sptr)
        .def("get_chef_lattice", &Lattice_simulator::get_chef_lattice_sptr)
        .def("calculate_lattice_functions",
                &Lattice_simulator::calculate_lattice_functions)
        ;

    class_<Step, Step_sptr >("Step", init<double >())
//            .def("append", remember how to overload methods...)
            .def("apply",&Step::apply)
            .def("get_operators",&Step::get_operators,
                    return_value_policy<copy_const_reference >())
            ;
    to_python_converter<Steps,
             container_conversions::to_tuple<Steps > >();

    class_<Stepper >("Stepper",no_init)
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

    void (Propagator::*propagate1)(Bunch &, int, Diagnostics &,
            Diagnostics &, bool) = &Propagator::propagate;
    void (Propagator::*propagate2)(Bunch &, int, Multi_diagnostics &,
            Multi_diagnostics &, bool) = &Propagator::propagate;

    class_<Propagator >("Propagator",init<Stepper_sptr >())
            .def("propagate", propagate1,
                    propagate_member_overloads())
            .def("propagate", propagate2,
                    propagate_member_overloads())
            ;
}
