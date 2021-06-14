#include "independent_operator.h"
#include "aperture_operation.h"
#include "synergia/utils/simple_timer.h"

void
Independent_operator::update_operations(
        Reference_particle const& reference_particle)
{
    operations.clear();
    operations_revisions.clear();

    std::string aperture_type("");
    bool need_left_aperture, need_right_aperture;
    std::string extractor_type(""), last_extractor_type("");
    // Group slices of equal extractor_type and pass to operation_extractor
    // to get operations.
    Lattice_element_slices group;
    for (Lattice_element_slices::iterator it = slices.begin();
            it != slices.end(); ++it) {
        if ((*it)->get_lattice_element().has_string_attribute(
                "aperture_type")) {
            aperture_type = (*it)->get_lattice_element().get_string_attribute(
                    "aperture_type");
            need_left_aperture = (*it)->has_left_edge();
            need_right_aperture = (*it)->has_right_edge();
        } else {
            need_left_aperture = false;
            need_right_aperture = false;
        }
        if ((*it)->get_lattice_element().has_string_attribute(
                "extractor_type")) {
            extractor_type = (*it)->get_lattice_element().get_string_attribute(
                    "extractor_type");
        } else {
            extractor_type = "default";
        }
        if (((extractor_type != last_extractor_type) || need_left_aperture)
                && (!group.empty())) {
            Independent_operations group_operations =
                    operation_extractor_map_sptr->get_extractor(last_extractor_type)->extract(
                            reference_particle, group);
            operations.splice(operations.end(), group_operations);
            group.clear();
        }
        if (need_left_aperture) {
            Aperture_operation_extractor_sptr extractor(
                    aperture_operation_extractor_map_sptr->get_extractor(
                            aperture_type));
            Aperture_operation_sptr aperture_operation_sptr(
                    extractor->extract(*it));
            operations.push_back(aperture_operation_sptr);

        }
        group.push_back(*it);
        last_extractor_type = extractor_type;
        if (need_right_aperture) {
            Independent_operations group_operations =
                    operation_extractor_map_sptr->get_extractor(extractor_type)->extract(
                            reference_particle, group);
            operations.splice(operations.end(), group_operations);
            group.clear();

            Aperture_operation_extractor_sptr extractor(
                    aperture_operation_extractor_map_sptr->get_extractor(
                            aperture_type));
            Aperture_operation_sptr aperture_operation_sptr(
                    extractor->extract(*it));
            operations.push_back(aperture_operation_sptr);
        }
        operations_revisions.push_back(
                (*it)->get_lattice_element().get_revision());
    }
    if (!group.empty()) {
        Independent_operations group_operations =
                operation_extractor_map_sptr->get_extractor(extractor_type)->extract(
                        reference_particle, group);
        operations.splice(operations.end(), group_operations);
    }
    Aperture_operation_extractor_sptr extractor(
            aperture_operation_extractor_map_sptr->get_extractor("default"));
    Aperture_operation_sptr aperture_operation_sptr(
            extractor->extract(slices.back()));
    operations.push_back(aperture_operation_sptr);
    Aperture_operation_sptr finite_aperture_operation_sptr(
            new Finite_aperture_operation(slices.back()));
    operations.push_back(finite_aperture_operation_sptr);

    have_operations = true;
    operations_reference_particle = reference_particle;
}

bool
Independent_operator::need_update(Reference_particle const& reference_particle,
        int verbosity, Logger & logger)
{
    const double reference_particle_tolerance = 1.0e-8;
    bool retval;
    if (have_operations) {
        retval = false;
        if (reference_particle.equal(operations_reference_particle,
                reference_particle_tolerance)) {
            std::list<long int >::const_iterator rev_it =
                    operations_revisions.begin();
            for (Lattice_element_slices::const_iterator it = slices.begin();
                    it != slices.end(); ++it) {
                long int cached_revision = (*rev_it);
                long int revision = (*it)->get_lattice_element().get_revision();
                if (revision != cached_revision) {
                    if (verbosity > 4) {
                        logger
                                << "Independent_operator: needs update because lattice element "
                                << (*it)->get_lattice_element().get_name()
                                << " has changed" << std::endl;
                    }
                    retval = true;
                    break;
                }
                ++rev_it;
            }
        } else {
            if (verbosity > 4) {
                logger
                        << "Independent_operator: needs update because reference particle has changed"
                        << std::endl;
            }
            retval = true;
        }
    } else {
        if (verbosity > 4) {
            logger
                    << "Independent_operator: needs update because does not have operations"
                    << std::endl;
        }
        retval = true;
    }
    return retval;
}

const char Independent_operator::type_name[] = "independent";

Independent_operator::Independent_operator(std::string const& name,
        Operation_extractor_map_sptr operation_extractor_map_sptr,
        Aperture_operation_extractor_map_sptr aperture_operation_extractor_map_sptr) :
        Operator(name, type_name), operation_extractor_map_sptr(
                operation_extractor_map_sptr), aperture_operation_extractor_map_sptr(
                aperture_operation_extractor_map_sptr), have_operations(false)
{
}

Independent_operator::Independent_operator()
{
}

void
Independent_operator::append_slice(Lattice_element_slice_sptr slice_sptr)
{
    slices.push_back(slice_sptr);
    have_operations = false;
}

Lattice_element_slices const&
Independent_operator::get_slices() const
{
    return slices;
}

Independent_operations const&
Independent_operator::get_operations() const
{
    return operations;
}

Independent_operations &
Independent_operator::get_operations()
{
    return operations;
}

void
Independent_operator::apply(Bunch & bunch, double time_step, Step & step,
        int verbosity, Diagnosticss const& per_operation_diagnosticss,
        Logger & logger)
{
    double t_total = simple_timer_current();
    if (verbosity > 4) {
        logger << "Independent_operator: slice(s) = " << std::endl;
        for (Lattice_element_slices::const_iterator it = slices.begin();
                it != slices.end(); ++it) {
            logger << "    " << (*it)->as_string() << std::endl;
        }
    }
    double t = simple_timer_current();
    bool do_update = need_update(bunch.get_reference_particle(), verbosity,
            logger);
    t = simple_timer_show(t, "independent_operator_apply-test_update");
    if (do_update) {
        if (verbosity > 3) {
            logger << "Independent_operator: updating operations" << std::endl;
        }
        update_operations(bunch.get_reference_particle());
        t = simple_timer_show(t, "independent_operator_apply-update_operations");
    }
    for (Independent_operations::iterator it = operations.begin();
            it != operations.end(); ++it) {
        if (verbosity > 3) {
            logger << "Independent_operator: operation type = "
                    << (*it)->get_type() << std::endl;
        }
        (*it)->apply(bunch, verbosity, logger);
        std::string label(
                "independent_operator_apply-" + (*it)->get_type()
                        + "_operation_apply");
        t = simple_timer_show(t, label.c_str());
        for (Diagnosticss::const_iterator itd =
                per_operation_diagnosticss.begin();
                itd != per_operation_diagnosticss.end(); ++it) {
            (*itd)->update_and_write();
        }
        t = simple_timer_show(t, "diagnostics-operation");
    }
    bunch.update_total_num();
    t = simple_timer_show(t, "independent_operator_apply-bunch_update_total_num");
    t_total = simple_timer_show(t_total, "independent_operator_apply-total");
}

void
Independent_operator::apply(Bunch & bunch, double time_step, Step & step,
        int verbosity, Logger & logger, Multi_diagnostics & diagnostics)
{

    if (need_update(bunch.get_reference_particle(), verbosity, logger)) {
        update_operations(bunch.get_reference_particle());
    }

    double t;
    t = simple_timer_current();
    for (Independent_operations::iterator it = operations.begin();
            it != operations.end(); ++it) {
        //  std::cout<<" opertor.cc operator name="<<(*it)->get_type()<<std::endl;
        for (Multi_diagnostics::iterator itd = diagnostics.begin();
                itd != diagnostics.end(); ++itd) {

            (*itd)->update_and_write();
        }
        (*it)->apply(bunch, verbosity, logger);
    }
    t = simple_timer_show(t, "independent_operator_apply2");
}

void
Independent_operator::print() const
{
    Operator::print();
    int count = 0;
    for (Lattice_element_slices::const_iterator it = slices.begin();
            it != slices.end(); ++it) {
        ++count;
        std::cout << "slice " << count << ": ";
        (*it)->print();
    }
}

template<class Archive>
    void
    Independent_operator::serialize(Archive & ar, const unsigned int version)
    {
        ar &
        BOOST_SERIALIZATION_BASE_OBJECT_NVP(Operator);
        ar & BOOST_SERIALIZATION_NVP(slices);
        ar & BOOST_SERIALIZATION_NVP(operations);
        ar & BOOST_SERIALIZATION_NVP(operations_revisions);
        ar & BOOST_SERIALIZATION_NVP(operations_reference_particle);
        ar & BOOST_SERIALIZATION_NVP(operation_extractor_map_sptr);
        ar & BOOST_SERIALIZATION_NVP(aperture_operation_extractor_map_sptr);
        ar & BOOST_SERIALIZATION_NVP(have_operations);
    }

template
void
Independent_operator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Independent_operator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Independent_operator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Independent_operator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Independent_operator::~Independent_operator()
{

}
BOOST_CLASS_EXPORT_IMPLEMENT(Independent_operator)