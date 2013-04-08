#include "operator.h"
#include "lattice_simulator.h"
#include "populate_stationary.h"
#include "stepper.h"
#include "independent_stepper.h"
#include "independent_stepper_elements.h"
#include "split_operator_stepper.h"
#include "split_operator_stepper_elements.h"
#include "split_operator_stepper_choice.h"
#include "propagator.h"
#include "propagate_actions.h"
#include "diagnostics_actions.h"
#include "dense_mapping.h"
#include "resume.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include "synergia/utils/numpy_multi_ref_converter.h"


using namespace boost::python;

struct Propagate_actions_callback : Propagate_actions
{
    Propagate_actions_callback(PyObject *p) :
        Propagate_actions(), self(handle< > (borrowed(p)))
    {
    }
    Propagate_actions_callback(PyObject *p, const Propagate_actions& x) :
        Propagate_actions(x), self(handle< > (borrowed(p)))
    {
    }
    // default constructor needed for serialization
    Propagate_actions_callback(): Propagate_actions(), self()
    {
    }
    void
    first_action(Stepper & stepper, Bunch & bunch)
    {
        call_method<void > (self.ptr(), "first_action", boost::ref(stepper),
                boost::ref(bunch));
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
        call_method<void > (self.ptr(), "turn_end_action", boost::ref(stepper),
                boost::ref(bunch), turn_num);
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
        call_method<void > (self.ptr(), "step_end_action", boost::ref(stepper),
                boost::ref(step), boost::ref(bunch),
                turn_num, step_num);
    }
    static void
    default_step_end_action(Propagate_actions& self_, Stepper & stepper,
            Step & step, Bunch & bunch, int turn_num, int step_num)
    {
        self_.Propagate_actions::step_end_action(stepper, step, bunch,
                turn_num, step_num);
    }
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const
        {
            ar &
            BOOST_SERIALIZATION_BASE_OBJECT_NVP(Propagate_actions);
            std::string pickled_object(
                    extract<std::string > (
                            import("cPickle").attr("dumps")(self)));
            ar & BOOST_SERIALIZATION_NVP(pickled_object);
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
        {
            ar &
            BOOST_SERIALIZATION_BASE_OBJECT_NVP(Propagate_actions);
            std::string pickled_object;
            ar & BOOST_SERIALIZATION_NVP(pickled_object);
            str pickle_str(pickled_object);
            self = import("cPickle").attr("loads")(pickle_str);
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
private:
    object self;
};
BOOST_CLASS_EXPORT_GUID(Propagate_actions_callback, "Propagate_actions_callback")

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads24,
        Propagator::propagate, 2, 4);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads35,
        Propagator::propagate, 3, 5);

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(propagate_member_overloads45,
//        Propagator::propagate, 4, 5);

Independent_operator_sptr
as_independent_operator(Operator_sptr & operator_sptr)
{
    return boost::dynamic_pointer_cast<Independent_operator >(operator_sptr);
}

Fast_mapping_operation_sptr
as_fast_mapping_operation(Independent_operation_sptr & independent_operation_sptr)
{
    return boost::dynamic_pointer_cast<Fast_mapping_operation >(independent_operation_sptr);
}

void (Diagnostics_actions::*add_per_turn1)(Diagnostics_sptr, int)
                            = &Diagnostics_actions::add_per_turn;
void (Diagnostics_actions::*add_per_turn2)(Diagnostics_sptr,
        std::list<int > const&)
                            = &Diagnostics_actions::add_per_turn;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_per_turn_member_overloads12,
        Diagnostics_actions::add_per_turn, 1, 2)
void (Diagnostics_actions::*add_per_step1)(Diagnostics_sptr, int)
                            = &Diagnostics_actions::add_per_step;
void (Diagnostics_actions::*add_per_step2)(Diagnostics_sptr,
        std::list<int > const&, int)
                            = &Diagnostics_actions::add_per_step;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_per_step_member_overloads12,
        Diagnostics_actions::add_per_step, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_per_step_member_overloads23,
        Diagnostics_actions::add_per_step, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_per_forced_diagnostics_step_member_overloads12,
        Diagnostics_actions::add_per_forced_diagnostics_step, 1, 2)

void (Bunch_simulator::*bs_add_per_turn1)(Diagnostics_sptr, int)
                            = &Bunch_simulator::add_per_turn;
void (Bunch_simulator::*bs_add_per_turn2)(Diagnostics_sptr,
        std::list<int > const&)
                            = &Bunch_simulator::add_per_turn;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bs_add_per_turn_member_overloads12,
        Bunch_simulator::add_per_turn, 1, 2)
void (Bunch_simulator::*bs_add_per_step1)(Diagnostics_sptr, int)
                            = &Bunch_simulator::add_per_step;
void (Bunch_simulator::*bs_add_per_step2)(Diagnostics_sptr,
        std::list<int > const&, int)
                            = &Bunch_simulator::add_per_step;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bs_add_per_step_member_overloads12,
        Bunch_simulator::add_per_step, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bs_add_per_step_member_overloads23,
        Bunch_simulator::add_per_step, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bs_add_per_forced_diagnostics_step_member_overloads12,
                Bunch_simulator::add_per_forced_diagnostics_step, 1, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_horizontal_tune_overloads01,
        Lattice_simulator::get_horizontal_tune, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_vertical_tune_overloads01,
        Lattice_simulator::get_vertical_tune, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_horizontal_chromaticity_overloads01, 
	Lattice_simulator::get_horizontal_chromaticity, 0,1)				
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_vertical_chromaticity_overloads01, 
	Lattice_simulator::get_vertical_chromaticity, 0,1)  
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_slip_factor_overloads01,
	Lattice_simulator::get_slip_factor,0,1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_momentum_compaction_overloads01,	
	Lattice_simulator::get_momentum_compaction,0,1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_both_tunes_overloads01,
        Lattice_simulator::get_both_tunes, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(adjust_chromaticities_overloads46,
		 Lattice_simulator::adjust_chromaticities, 4, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(adjust_tunes_overloads45,
		 Lattice_simulator::adjust_tunes, 4, 5)



void (Bunch_train_simulator::*bts_add_per_turn1)(int, Diagnostics_sptr, int)
                            = &Bunch_train_simulator::add_per_turn;
			    
void (Bunch_train_simulator::*bts_add_per_turn2)(int, Diagnostics_sptr,
        std::list<int > const&)
                            = &Bunch_train_simulator::add_per_turn;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bts_add_per_turn_member_overloads12,
        Bunch_train_simulator::add_per_turn, 2, 3)
void (Bunch_train_simulator::*bts_add_per_step1)(int, Diagnostics_sptr, int)
                            = &Bunch_train_simulator::add_per_step;
void (Bunch_train_simulator::*bts_add_per_step2)(int, Diagnostics_sptr,
        std::list<int > const&, int)
                            = &Bunch_train_simulator::add_per_step;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bts_add_per_step_member_overloads12,
        Bunch_train_simulator::add_per_step, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bts_add_per_step_member_overloads23,
        Bunch_train_simulator::add_per_step, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(bts_add_per_forced_diagnostics_step_member_overloads12,
                Bunch_train_simulator::add_per_forced_diagnostics_step, 2, 3)

void
(Operator::*apply_bunch)(Bunch &, double, Step &, int, Diagnosticss const&,
        Logger &) = &Operator::apply;
void
(Operator::*apply_bunch_train)(Bunch_train &, double, Step &, int,
        Train_diagnosticss const&, Logger &) = &Operator::apply;

BOOST_PYTHON_MODULE(simulation)
{
    import_array();

    class_<Operator, Operator_sptr, boost::noncopyable >("Operator", no_init)
        .def("get_name", &Operator::get_name,
                return_value_policy<copy_const_reference >())
        .def("get_type", &Operator::get_type,
                return_value_policy<copy_const_reference >())
        .def("apply", apply_bunch)
        .def("apply", apply_bunch_train)
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

    container_conversions::from_python_sequence<Collective_operators,
             container_conversions::variable_capacity_policy >();

    class_<Dummy_collective_operator, Dummy_collective_operator_sptr,
        bases<Collective_operator> >("Dummy_collective_operator",
                init<std::string const& >())
//        .def("get_name", &Collective_operator::get_name)
//        .def("get_type", &Collective_operator::get_type)
//        .def("apply", &Collective_operator::apply)
//        .def("print_", &Collective_operator::print)
        ;

    Independent_operations &
    (Independent_operator::*get_operations1)() =
            &Independent_operator::get_operations;

    class_<Independent_operator, Independent_operator_sptr,
        bases<Operator > >("Independent_operator", init<std::string const&,
                Operation_extractor_map_sptr,
                Aperture_operation_extractor_map_sptr >())
//        .def("get_name", &Collective_operator::get_name)
//        .def("get_type", &Collective_operator::get_type)
//        .def("apply", &Collective_operator::apply)
//        .def("print_", &Collective_operator::print)
        .def("append_slice", &Independent_operator::append_slice)
        .def("get_slices", &Independent_operator::get_slices,
                return_value_policy<copy_const_reference >())
        .def("get_operations", get_operations1,
                return_value_policy<copy_non_const_reference >())
        ;

    def("as_independent_operator", as_independent_operator);

    to_python_converter<Operators,
             container_conversions::to_tuple<Operators > >();

    class_<Independent_operation, Independent_operation_sptr,
        boost::noncopyable >("Independent_operation", no_init)
        .def("get_type", &Independent_operation::get_type,
                return_value_policy<copy_const_reference >())
        .def("apply", &Independent_operation::apply)
        ;

    to_python_converter<Independent_operations,
             container_conversions::to_tuple<Independent_operations > >();

    class_<Fast_mapping >("Fast_mapping", init<int >())
            .def(init<std::string const& >())
            .def("set_length", &Fast_mapping::set_length)
            .def("get_length", &Fast_mapping::get_length)
            .def("apply", &Fast_mapping::apply)
            .def("as_string", &Fast_mapping::as_string)
            .def("write_to_file", &Fast_mapping::write_to_file)
            ;

    class_<Dense_mapping >("Dense_mapping", init<Fast_mapping const& >())
            .def("get_constant_term", &Dense_mapping::get_constant_term)
            .def("get_linear_term", &Dense_mapping::get_linear_term)
            ;

    class_<Fast_mapping_operation, Fast_mapping_operation_sptr,
        bases<Independent_operation > >("Fast_mapping_operation", no_init)
        .def("get_fast_mapping", &Fast_mapping_operation::get_fast_mapping,
                return_value_policy<copy_const_reference >())
        ;

    def("as_fast_mapping_operation", as_fast_mapping_operation);

    class_<Chef_propagate_operation, Chef_propagate_operation_sptr,
        bases<Independent_operation > >("Chef_propagate_operation", no_init)
        ;

    Lattice_functions const&
    (Lattice_simulator::*get_lattice_functions1)(Lattice_element &) =
            &Lattice_simulator::get_lattice_functions;
    Lattice_functions const&
    (Lattice_simulator::*get_lattice_functions2)(Lattice_element_slice &) =
            &Lattice_simulator::get_lattice_functions;

    class_<Aperture_operation_extractor, boost::noncopyable >
            Aperture_operation_extractor_("Aperture_operation_extractor_",
                    no_init);
    Aperture_operation_extractor_.def("extract",
            &Aperture_operation_extractor::extract);

    class_<Circular_extractor, bases<Aperture_operation_extractor > >
            Circular_extractor_("Circular_extractor", init< > ());
    Circular_extractor_.def("extract", &Circular_extractor::extract);

    class_<Elliptical_extractor, bases<Aperture_operation_extractor > >
            Elliptical_extractor_("Elliptical_extractor", init< > ());
    Elliptical_extractor_.def("extract", &Elliptical_extractor::extract);

    class_<Rectangular_extractor, bases<Aperture_operation_extractor > >
            Rectangular_extractor_("Rectangular_extractor", init< > ());
    Rectangular_extractor_.def("extract", &Rectangular_extractor::extract);

    class_<Polygon_extractor, bases<Aperture_operation_extractor > >
            Polygon_extractor_("Polygon_extractor", init< > ());
    Polygon_extractor_.def("extract", &Polygon_extractor::extract);

    class_<Wire_elliptical_extractor, bases<Aperture_operation_extractor > >
            Wire_elliptical_extractor_("Wire_elliptical_extractor", init< > ());
    Wire_elliptical_extractor_.def("extract", &Wire_elliptical_extractor::extract);

    class_<Lambertson_extractor, bases<Aperture_operation_extractor > >
            Lambertson_extractor_("Lambertson_extractor", init< > ());
    Lambertson_extractor_.def("extract", &Lambertson_extractor::extract);

    class_<Aperture_operation_extractor_map,
            Aperture_operation_extractor_map_sptr >
            Aperture_operation_extractor_map_(
                    "Aperture_operation_extractor_map", init< > ());
    Aperture_operation_extractor_map_.def("set_extractor",
            &Aperture_operation_extractor_map::set_extractor);
    Aperture_operation_extractor_map_.def("get_extractor",
            &Aperture_operation_extractor_map::get_extractor);
    Aperture_operation_extractor_map_.def("get_extractor_names",
            &Aperture_operation_extractor_map::get_extractor_names);

    to_python_converter <std::pair<double,double >,
    container_conversions::PairToTupleConverter<double, double > >();

    class_<Lattice_simulator >("Lattice_simulator",
            init<Lattice_sptr, int >())
        .def("set_slices",
                &Lattice_simulator::set_slices)
        .def("get_slices",
                &Lattice_simulator::get_slices,
                return_value_policy<copy_const_reference >())
        .def("get_map_order", &Lattice_simulator::get_map_order)
        .def("get_operation_extractor_map",
                &Lattice_simulator::get_operation_extractor_map_sptr)
        .def("get_aperture_operation_extractor_map",
                &Lattice_simulator::get_aperture_operation_extractor_map_sptr)
        .def("get_lattice", &Lattice_simulator::get_lattice_sptr)
        .def("get_chef_lattice", &Lattice_simulator::get_chef_lattice_sptr)
        .def("get_bucket_length", &Lattice_simulator::get_bucket_length)
        .def("get_number_buckets",&Lattice_simulator::get_number_buckets)
        .def("update", &Lattice_simulator::update)
        .def("calculate_element_lattice_functions",
                &Lattice_simulator::calculate_element_lattice_functions)
        .def("calculate_slice_lattice_functions",
                &Lattice_simulator::calculate_slice_lattice_functions)
        .def("get_lattice_functions", get_lattice_functions1,
                return_value_policy<copy_const_reference >())
        .def("get_lattice_functions", get_lattice_functions2,
                return_value_policy<copy_const_reference >())
        .def("get_horizontal_tune", &Lattice_simulator::get_horizontal_tune,
	                          get_horizontal_tune_overloads01())
        .def("get_vertical_tune", &Lattice_simulator::get_vertical_tune,
	                          get_vertical_tune_overloads01())
        .def("get_both_tunes", &Lattice_simulator::get_both_tunes,
	                          get_both_tunes_overloads01())
        .def("adjust_tunes", &Lattice_simulator::adjust_tunes,
	                           adjust_tunes_overloads45())
        .def("get_horizontal_chromaticity", &Lattice_simulator::get_horizontal_chromaticity,
	                     get_horizontal_chromaticity_overloads01())
        .def("get_vertical_chromaticity", &Lattice_simulator::get_vertical_chromaticity,
	                     get_vertical_chromaticity_overloads01())
        .def("get_momentum_compaction", &Lattice_simulator::get_momentum_compaction,
                                  get_momentum_compaction_overloads01())
        .def("get_slip_factor",&Lattice_simulator::get_slip_factor,
	                               get_slip_factor_overloads01())
        .def("adjust_chromaticities", &Lattice_simulator::adjust_chromaticities,
	                             adjust_chromaticities_overloads46())
      .def("is_ring", &Lattice_simulator::is_ring)
      .def("get_linear_one_turn_map", &Lattice_simulator::get_linear_one_turn_map)
      .def("check_linear_normal_form", &Lattice_simulator::check_linear_normal_form)
      .def("convert_normal_to_human", &Lattice_simulator::convert_normal_to_human)
      .def("convert_human_to_normal", &Lattice_simulator::convert_human_to_normal)
      .def("get_stationary_actions", &Lattice_simulator::get_stationary_actions)
//         .def("print_CS_lattice_functions", &Lattice_simulator::print_cs_lattice_functions)
//         .def("print_ET_lattice_functions", &Lattice_simulator::print_et_lattice_functions)
//         .def("print_LB_lattice_functions", &Lattice_simulator::print_lb_lattice_functions)
//         .def("print_dispersion_closedOrbit", &Lattice_simulator::print_dispersion_closedOrbit)
        .def("print_lattice_functions", &Lattice_simulator::print_lattice_functions)
        ;

    def("populate_6d_stationary_torus", &populate_6d_stationary_torus);
    def("populate_6d_stationary_gaussian", &populate_6d_stationary_gaussian);
    def("populate_6d_stationary_truncated_longitudinal_gaussian", &populate_6d_stationary_truncated_longitudinal_gaussian);
    def("populate_6d_stationary_clipped_longitudinal_gaussian", &populate_6d_stationary_clipped_longitudinal_gaussian);

    class_<Lattice_functions >("Lattice_functions",
            init<LattFuncSage::lattFunc const& >())
        .def_readonly("alpha_x", &Lattice_functions::alpha_x)
        .def_readonly("alpha_y", &Lattice_functions::alpha_y)
        .def_readonly("beta_x", &Lattice_functions::beta_x)
        .def_readonly("beta_y", &Lattice_functions::beta_y)
        .def_readonly("psi_x", &Lattice_functions::psi_x)
        .def_readonly("psi_y", &Lattice_functions::psi_y)
        .def_readonly("D_x", &Lattice_functions::D_x)
        .def_readonly("D_y", &Lattice_functions::D_y)
        .def_readonly("Dprime_x", &Lattice_functions::Dprime_x)
        .def_readonly("Dprime_y", &Lattice_functions::Dprime_y)
        .def_readonly("arc_length", &Lattice_functions::arc_length)
        ;

    void (Step::*apply1)(Bunch &, int, Diagnosticss const&, Diagnosticss const&, Logger &) = &Step::apply;
    void (Step::*apply2)(Bunch_train &, int, Train_diagnosticss const&, Train_diagnosticss const&, Logger &) = &Step::apply;
    Operators& (Step::*get_operators1)() = &Step::get_operators;
    class_<Step, Step_sptr >("Step", init<double >())
//            .def("append", remember how to overload methods...)
            .def("apply",apply1)
            .def("apply",apply2)
            .def("get_operators",get_operators1,
                    return_value_policy<copy_non_const_reference >())
            .def("get_time_fractions",&Step::get_time_fractions,
                    return_value_policy<copy_const_reference >())
            .def("print_", &Step::print)
            .def("get_length", &Step::get_length)
            ;

    to_python_converter<Steps,
             container_conversions::to_tuple<Steps > >();

    class_<Stepper >("Stepper",init<Lattice_sptr, int >())
        .def(init<Lattice_simulator const& >())
        .def("get_lattice_simulator", &Stepper::get_lattice_simulator,
                return_internal_reference< >())
        .def("get_steps", &Stepper::get_steps,
                return_value_policy<copy_non_const_reference >())
        .def("force_update_operations_no_collective",
                &Stepper::force_update_operations_no_collective)
        .def("print_", &Stepper::print)
        ;

    class_<Split_operator_stepper, bases<Stepper > >("Split_operator_stepper",
            init<Lattice_sptr, int, Collective_operator_sptr, int >())
           .def(init<Lattice_sptr, int, Collective_operators, int >())
           .def(init<Lattice_simulator const&, Collective_operator_sptr, int >())
           .def(init<Lattice_simulator const&, Collective_operators, int >())
            ;

    class_<Split_operator_stepper_elements, bases<Stepper > >("Split_operator_stepper_elements",
            init<Lattice_sptr, int, Collective_operator_sptr, int >())
	    .def(init<Lattice_sptr, int, Collective_operators const &, int  >())
	    .def(init<Lattice_simulator const&, Collective_operator_sptr, int >())
	    .def(init<Lattice_simulator const&, Collective_operators const &, int  >())
	    ;

    class_<Kicks >("Kicks", init< >())
            .def(init<Collective_operators const& , int >())
            ;

    class_<List_choice_map >("List_choice_map")
        .def(map_indexing_suite<std::map<std::string, Kicks> >() );
        ;

    class_<Split_operator_stepper_choice, bases<Stepper > >("Split_operator_stepper_choice",
            init<Lattice_sptr, int, List_choice_map const&, optional<bool > >())
           .def(init< int, Lattice_sptr, int, List_choice_map const&, optional<bool > >())
           .def(init<Lattice_simulator const&, List_choice_map const&, optional<bool > >())
           .def(init< int, Lattice_simulator const&, List_choice_map const&, optional<bool > >())
           //.def(init<Lattice_simulator const&, List_choice_map const&, int>())
         //  .def(init<Lattice_simulator const&, List_choice_map const& >())
           ;

    class_<Independent_stepper, bases<Stepper > >("Independent_stepper",
            init<Lattice_sptr, int, int >())
            .def(init<Lattice_simulator const&, int >())
            ;

    class_<Independent_stepper_elements, bases<Stepper > >("Independent_stepper_elements",
            init<Lattice_sptr, int, int >())
            .def(init<Lattice_simulator const&, int >())
            ;

    class_<Propagate_actions, Propagate_actions_callback >("Propagate_actions",
            init< >())
            .def("first_action",
                    &Propagate_actions_callback::default_first_action)
            .def("turn_end_action",
                    &Propagate_actions_callback::default_turn_end_action)
            .def("step_end_action",
                    &Propagate_actions_callback::default_step_end_action)
            .enable_pickling()
            ;

    class_<Diagnostics_actions, Diagnostics_actions_sptr >(
            "Diagnostics_actions", init< >())
            .def("add_per_turn", add_per_turn1,
                    add_per_turn_member_overloads12())
            .def("add_per_turn", add_per_turn2)
            .def("add_per_step", add_per_step1,
                    add_per_turn_member_overloads12())
            .def("add_per_step", add_per_step2,
                    add_per_step_member_overloads23())
            .def("add_per_forced_diagnostics_step",
                    &Diagnostics_actions::add_per_forced_diagnostics_step,
                    add_per_forced_diagnostics_step_member_overloads12())
            .def("first_action", &Diagnostics_actions::first_action)
            .def("turn_end_action", &Diagnostics_actions::turn_end_action)
            .def("step_end_action", &Diagnostics_actions::step_end_action)
            ;

    class_<Bunch_simulator > ("Bunch_simulator", init<Bunch_sptr>())
            .def(init<Bunch_sptr, Diagnostics_actions_sptr>())
            .def("get_bunch", &Bunch_simulator::get_bunch_sptr)
            .def("get_diagnostics_actions", &Bunch_simulator::get_diagnostics_actions_sptr)
            .def("add_per_turn", bs_add_per_turn1,
                    bs_add_per_turn_member_overloads12())
            .def("add_per_turn", bs_add_per_turn2)
            .def("add_per_step", bs_add_per_step1,
                    bs_add_per_step_member_overloads12())
            .def("add_per_step", bs_add_per_step2,
                    bs_add_per_step_member_overloads23())
            .def("add_per_forced_diagnostics_step",
                    &Bunch_simulator::add_per_forced_diagnostics_step,
                    bs_add_per_forced_diagnostics_step_member_overloads12())
            ;

    class_<Bunch_train_simulator > ("Bunch_train_simulator", init<Bunch_train_sptr>())
            .def("get_bunch_train", &Bunch_train_simulator::get_bunch_train_sptr)
            .def("get_diagnostics_actionss", &Bunch_train_simulator::get_diagnostics_actionss,
                    return_internal_reference< >())
            .def("add_per_turn", bts_add_per_turn1,
                    bts_add_per_turn_member_overloads12())
            .def("add_per_turn", bts_add_per_turn2)
            .def("add_per_step", bts_add_per_step1,
                    bts_add_per_step_member_overloads12())
            .def("add_per_step", bts_add_per_step2,
                    bts_add_per_step_member_overloads23())
            .def("add_per_forced_diagnostics_step",
                    &Bunch_train_simulator::add_per_forced_diagnostics_step,
                    bts_add_per_forced_diagnostics_step_member_overloads12())
            ;

    void (Propagator::*propagate1)(Bunch_simulator &, int, int, int)
                                = &Propagator::propagate;

    void (Propagator::*propagate2)(Bunch_simulator &, Propagate_actions &, int, int, int)
                                = &Propagator::propagate;

    void (Propagator::*propagate_train1)(Bunch_train_simulator &, int, int, int)
                                = &Propagator::propagate;

    void (Propagator::*propagate_train2)(Bunch_train_simulator &, Propagate_actions &, int, int, int)
                                = &Propagator::propagate;

//    void (Propagator::*propagate3)(Bunch_with_diagnostics_train &, int, Propagate_actions &, bool)
//                                = &Propagator::propagate;
//
//    void (Propagator::*propagate4)(Bunch_with_diagnostics_train &, int, bool)
//                                = &Propagator::propagate;

//    void (Propagator::*propagate5)(Bunch &, int, Diagnostics &,
//            Diagnostics &, bool) = &Propagator::propagate;

//    void (Propagator::*propagate6)(Bunch &, int, Multi_diagnostics &,
//            Multi_diagnostics &, bool) = &Propagator::propagate;
//
//    void (Propagator::*propagate7)(Bunch &, int, Diagnostics_actions&,
//            int) = &Propagator::propagate;
//
//    void (Propagator::*propagate8)(Bunch &, int, Diagnostics_actions &,
//            Propagate_actions &, int) = &Propagator::propagate;

    {
    scope propagator_scope(
    class_<Propagator >("Propagator",init<Stepper_sptr >())
            .def_readonly("default_checkpoint_dir", &Propagator::default_checkpoint_dir)
            .def_readonly("description_file_name", &Propagator::description_file_name)
            .def_readonly("propagator_archive_name", &Propagator::propagator_archive_name)
            .def_readonly("propagator_xml_archive_name", &Propagator::propagator_xml_archive_name)
            .def_readonly("state_archive_name", &Propagator::state_archive_name)
            .def_readonly("state_xml_archive_name", &Propagator::state_xml_archive_name)
            .def_readonly("log_file_name", &Propagator::log_file_name)
            .def_readonly("stop_file_name", &Propagator::stop_file_name)
            .def_readonly("alt_stop_file_name", &Propagator::alt_stop_file_name)
            .def("get_stepper", &Propagator::get_stepper_sptr)
            .def("set_checkpoint_period", &Propagator::set_checkpoint_period)
            .def("get_checkpoint_period", &Propagator::get_checkpoint_period)
            .def("set_checkpoint_dir", &Propagator::set_checkpoint_dir)
            .def("get_checkpoint_dir", &Propagator::get_checkpoint_dir,
                    return_value_policy<copy_const_reference >())
            .def("set_checkpoint_with_xml", &Propagator::set_checkpoint_with_xml)
            .def("get_checkpoint_with_xml", &Propagator::get_checkpoint_with_xml)
            .def("set_final_checkpoint", &Propagator::set_final_checkpoint)
            .def("get_final_checkpoint", &Propagator::get_final_checkpoint)
            .def("set_concurrent_io", &Propagator::set_concurrent_io)
            .def("get_concurrent_io", &Propagator::get_concurrent_io)
            .def("propagate", propagate1,
                 propagate_member_overloads24())
            .def("propagate", propagate2,
                    propagate_member_overloads35())
            .def("propagate", propagate_train1,
                 propagate_member_overloads24())
            .def("propagate", propagate_train2,
                    propagate_member_overloads35())
            .def("get_resume_state", &Propagator::get_resume_state)
            .def("resume", &Propagator::resume)
//            .def("resume", &Propagator::resume,
//                    resume_member_overloads())
//            .def("propagate", propagate3,
//                    propagate_member_overloads34())
//            .def("propagate", propagate4,
//                    propagate_member_overloads23())
//            .def("propagate", propagate5,
//                    propagate_member_overloads45())
//             .def("propagate", propagate6,
//                    propagate_member_overloads45())
//            .def("propagate", propagate7,
//                    propagate_member_overloads34())
//            .def("propagate", propagate8,
//                    propagate_member_overloads45())
    );

    class_<Propagator::State>("State", no_init)
            .def_readwrite("bunch_simulator", &Propagator::State::bunch_simulator_ptr)
            .def_readwrite("propagate_actions", &Propagator::State::propagate_actions_ptr)
            .def_readwrite("num_turns", &Propagator::State::num_turns)
            .def_readwrite("first_turn", &Propagator::State::first_turn)
            .def_readwrite("max_turns", &Propagator::State::max_turns)
            .def_readwrite("verbosity", &Propagator::State::verbosity)
            ;

    }

    {
    scope resume_scope(
    class_<Resume >("Resume", init<>())
            .def(init<std::string const& >())
            .def("set_checkpoint_period", &Resume::set_checkpoint_period)
            .def("get_checkpoint_period", &Resume::get_checkpoint_period)
            .def("set_new_checkpoint_dir", &Resume::set_new_checkpoint_dir)
            .def("get_new_checkpoint_dir", &Resume::get_new_checkpoint_dir,
                    return_value_policy<copy_const_reference >())
            .def("get_content", &Resume::get_content)
            .def("propagate", &Resume::propagate)
    );

    class_<Resume::Content>("Content", no_init)
            .def_readwrite("bunch", &Resume::Content::bunch_sptr)
            .def_readwrite("stepper", &Resume::Content::stepper_sptr)
            .def_readwrite("lattice", &Resume::Content::lattice_sptr)
            ;
    }
}
