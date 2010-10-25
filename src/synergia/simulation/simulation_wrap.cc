#include "operator.h"
#include "lattice_simulator.h"
#include "stepper.h"
#include "propagator.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

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
        bases<Operator> >("Collective_operator", init<std::string const& >())
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
        ;

//    class Step
//    {
//    private:
//        Operators operators;
//    public:
//        Step();
//        void
//        append(Operator_sptr operator_sptr);
//        void
//        append(Operators const& operators);
//        void
//        apply(Bunch & bunch);
//        Operators const&
//        get_operators() const;
//        void
//        print(int index) const;
//    };
//
//    typedef boost::shared_ptr<Step > Step_sptr;
//    typedef std::list<Step_sptr > Steps;
    class_<Step, Step_sptr >("Step", init<>())
//            .def("append", remember how to overload methods...)
            .def("apply",&Step::apply)
            .def("get_operators",&Step::get_operators,
                    return_value_policy<copy_const_reference >())
            ;
    to_python_converter<Steps,
             container_conversions::to_tuple<Steps > >();

    //    class Stepper
//    {
//    private:
//        Steps steps;
//
//    public:
//        Steps &
//        get_steps();
//        virtual void
//        print() const;
//
//        virtual
//        ~Stepper();
//    };

    class_<Stepper >("Stepper",no_init)
        .def("get_steps", &Stepper::get_steps,
                return_value_policy<copy_non_const_reference >())
        .def("print_", &Stepper::print)
        ;

//    class Split_operator_stepper : public Stepper
//    {
//    private:
//        Lattice_simulator lattice_simulator;
//        Independent_operator_sptr
//        get_half_step(std::string const& name,
//                Lattice_elements::iterator & lattice_it, double & left,
//                Lattice_elements::iterator const & lattice_end,
//                const double half_step_length);
//        void
//        construct(Collective_operators const & collective_operators, int num_steps);
//    public:
//        Split_operator_stepper(Lattice_simulator const& lattice_simulator,
//                Collective_operator_sptr collective_operator, int num_steps);
//        Split_operator_stepper(Lattice_simulator const& lattice_simulator,
//                Collective_operators const & collective_operators, int num_steps);
//        ~Split_operator_stepper();
//    };
    class_<Split_operator_stepper, bases<Stepper > >("Split_operator_stepper",
            init<Lattice_simulator const&, Collective_operator_sptr, int >());

//    class Propagator
//    {
//    private:
//        Stepper_sptr stepper_sptr;
//
//        void
//        construct();
//    public:
//        Propagator(Stepper_sptr stepper_sptr);
//        void
//        propagate(Bunch & bunch, int num_turns, bool diagnostics_per_step,
//                bool diagnostics_per_turn);
//        ~Propagator();
//    };
    class_<Propagator >("Propagator",init<Stepper_sptr >())
            .def("propagate",&Propagator::propagate)
            ;
}
