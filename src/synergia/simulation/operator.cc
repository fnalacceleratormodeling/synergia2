#include <iostream>

#include "operator.h"
#include "synergia/simulation/operation_extractor.h"
#include "synergia/simulation/aperture_operation.h"



#if 0
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
    for (auto const & slice : slices)
    {
        auto const & element = slice.get_lattice_element();

        if (element.has_string_attribute("aperture_type")) 
        {
            aperture_type = element.get_string_attribute("aperture_type");
            need_left_aperture = slice.has_left_edge();
            need_right_aperture = slice.has_right_edge();
        } 
        else 
        {
            need_left_aperture = false;
            need_right_aperture = false;
        }

        extractor_type = element.get_string_attribute("extractor_type", "default");

        if (((extractor_type != last_extractor_type) || need_left_aperture) 
                && (!group.empty())) 
        {
            auto group_operations = extract_independent_operations(
                    extractor_type,
                    chef_lattice,
                    reference_particle,
                    group,
                    map_order);

            operations.push_back(group_operations);
            group.clear();
        }

        if (need_left_aperture) 
        {
            auto aperture_opn = extract_aperture_operation(aperture_type, slice);
            operations.push_back(aperture_opn);
        }

        group.push_back(slice);
        last_extractor_type = extractor_type;

        if (need_right_aperture) 
        {
            auto group_operations = extract_independent_operations(
                    extractor_type, chef_lattice,
                    reference_particle, group, map_order);

            operations.push_back(group_operations);
            group.clear();

            auto aperture_opn = extract_aperture_operation(aperture_type, slice);
            operations.push_back(aperture_opn);
        }

        operations_revisions.push_back(element.get_revision());
    }

    if (!group.empty()) 
    {
        auto group_operations = extract_independent_operations(
                extractor_type, chef_lattice,
                reference_particle, group, map_order);

        operations.push_back(group_operations);
    }

    auto aperture_opn = extract_aperture_operation("default", slices.back());
    operations.push_back(aperture_opn);

    auto finite_aperture = std::make_unique<Finite_aperture_operation>(slices.back());
    operations.push_back(finite_aperture);

    have_operations = true;
    operations_reference_particle = reference_particle;
}

bool Independent_operator::need_update(
        Reference_particle const& reference_particle,
        int verbosity, Logger & logger)
{
    const double reference_particle_tolerance = 1.0e-8;
    bool retval;

    if (have_operations) 
    {
        retval = false;

        if (reference_particle.equal(operations_reference_particle,
                reference_particle_tolerance)) 
        {
            std::list<long int >::const_iterator rev_it =
                    operations_revisions.begin();

            for (Lattice_element_slices::const_iterator it = slices.begin();
                    it != slices.end(); ++it) 
            {
                long int cached_revision = (*rev_it);
                long int revision = (*it)->get_lattice_element().get_revision();

                if (revision != cached_revision) 
                {
                    if (verbosity > 4) {
                        logger << "Independent_operator: needs update because lattice element "
                               << (*it)->get_lattice_element().get_name()
                               << " has changed" << std::endl;
                    }
                    retval = true;
                    break;
                }

                ++rev_it;
            }

        } 
        else 
        {
            if (verbosity > 4) {
                logger << "Independent_operator: needs update because reference particle has changed"
                       << std::endl;
            }

            retval = true;
        }

    } 
    else 
    {
        if (verbosity > 4) {
            logger << "Independent_operator: needs update because does not have operations"
                   << std::endl;
        }
        retval = true;
    }

    return retval;
}
#endif

Independent_operator::Independent_operator(std::string const & name, double time)
    : Operator(name, "independent", time)
    , slices()
    , operations()
{
}

bool Independent_operator::need_update(
        Reference_particle const & ref,
        Logger & logger)
{
    const double reference_particle_tolerance = 1.0e-8;
    bool retval = false;

#if 0
    if (have_operations) 
    {
        retval = false;

        if (reference_particle.equal(operations_reference_particle,
                reference_particle_tolerance)) 
        {
            std::list<long int >::const_iterator rev_it =
                    operations_revisions.begin();

            for (Lattice_element_slices::const_iterator it = slices.begin();
                    it != slices.end(); ++it) 
            {
                long int cached_revision = (*rev_it);
                long int revision = (*it)->get_lattice_element().get_revision();

                if (revision != cached_revision) 
                {
                    if (verbosity > 4) {
                        logger << "Independent_operator: needs update because lattice element "
                               << (*it)->get_lattice_element().get_name()
                               << " has changed" << std::endl;
                    }
                    retval = true;
                    break;
                }

                ++rev_it;
            }

        } 
        else 
        {
            if (verbosity > 4) {
                logger << "Independent_operator: needs update because reference particle has changed"
                       << std::endl;
            }

            retval = true;
        }

    } 
    else 
    {
        if (verbosity > 4) {
            logger << "Independent_operator: needs update because does not have operations"
                   << std::endl;
        }
        retval = true;
    }
#endif

    return retval;
}

void
Independent_operator::update_operations(
        Reference_particle const & reference_particle)
{
#if 0
    operations.clear();
    operations_revisions.clear();

    std::string aperture_type("");
    bool need_left_aperture, need_right_aperture;

    std::string extractor_type(""), last_extractor_type("");

    // Group slices of equal extractor_type and pass to operation_extractor
    // to get operations.
    Lattice_element_slices group;
    for (auto const & slice : slices)
    {
        auto const & element = slice.get_lattice_element();

        if (element.has_string_attribute("aperture_type")) 
        {
            aperture_type = element.get_string_attribute("aperture_type");
            need_left_aperture = slice.has_left_edge();
            need_right_aperture = slice.has_right_edge();
        } 
        else 
        {
            need_left_aperture = false;
            need_right_aperture = false;
        }

        extractor_type = element.get_string_attribute("extractor_type", "default");

        if (((extractor_type != last_extractor_type) || need_left_aperture) 
                && (!group.empty())) 
        {
            auto group_operations = extract_independent_operations(
                    extractor_type,
                    chef_lattice,
                    reference_particle,
                    group,
                    map_order);

            operations.push_back(group_operations);
            group.clear();
        }

        if (need_left_aperture) 
        {
            auto aperture_opn = extract_aperture_operation(aperture_type, slice);
            operations.push_back(aperture_opn);
        }

        group.push_back(slice);
        last_extractor_type = extractor_type;

        if (need_right_aperture) 
        {
            auto group_operations = extract_independent_operations(
                    extractor_type, chef_lattice,
                    reference_particle, group, map_order);

            operations.push_back(group_operations);
            group.clear();

            auto aperture_opn = extract_aperture_operation(aperture_type, slice);
            operations.push_back(aperture_opn);
        }

        operations_revisions.push_back(element.get_revision());
    }

    if (!group.empty()) 
    {
        auto group_operations = extract_independent_operations(
                extractor_type, chef_lattice,
                reference_particle, group, map_order);

        operations.push_back(group_operations);
    }

    auto aperture_opn = extract_aperture_operation("default", slices.back());
    operations.push_back(aperture_opn);

    auto finite_aperture = std::make_unique<Finite_aperture_operation>(slices.back());
    operations.push_back(finite_aperture);

    have_operations = true;
    operations_reference_particle = reference_particle;
#endif
}

void
Independent_operator::create_operations_impl(
        Lattice const & lattice)
{
    operations.clear();
    //operations_revisions.clear();

    std::string aperture_type("");
    bool need_left_aperture, need_right_aperture;

    std::string extractor_type(""), last_extractor_type("");

    // Group slices of equal extractor_type and pass to operation_extractor
    // to get operations.
    std::vector<Lattice_element_slice> group;
    for (auto const & slice : slices)
    {
        auto const & element = slice.get_lattice_element();

        if (element.has_string_attribute("aperture_type")) 
        {
            aperture_type = element.get_string_attribute("aperture_type");
            need_left_aperture = slice.has_left_edge();
            need_right_aperture = slice.has_right_edge();
        } 
        else 
        {
            need_left_aperture = false;
            need_right_aperture = false;
        }

        extractor_type = element.get_string_attribute("extractor_type", "default");

        if ( ((extractor_type != last_extractor_type) || need_left_aperture) 
                && (!group.empty()) ) 
        {
            extract_independent_operations(extractor_type, lattice, group, operations);
            group.clear();
        }

        if (need_left_aperture) 
        {
            operations.emplace_back(extract_aperture_operation(aperture_type, slice));
        }

        group.push_back(slice);
        last_extractor_type = extractor_type;

        if (need_right_aperture) 
        {
            extract_independent_operations(extractor_type, lattice, group, operations);
            operations.emplace_back(extract_aperture_operation(aperture_type, slice));
            group.clear();
        }

        //operations_revisions.push_back(element.get_revision());
    }

    if (!group.empty()) 
    {
        extract_independent_operations(extractor_type, lattice, group, operations);
    }

    // only the finite aperture will be attached to an independent operator by default
    operations.emplace_back(
            extract_aperture_operation(Finite_aperture::type, slices.back()));

#if 0
    have_operations = true;
    operations_reference_particle = reference_particle;
#endif
}


void
Independent_operator::apply_impl(
        Bunch_simulator & simulator, 
        double time_step,
        Logger & logger)
{
    using LV = LoggerV;

    //double t_total = simple_timer_current();

    logger(LV::DINFO) << "    Independent_operator: slice(s) = ";

    for (auto const & slice : slices)
        logger(LV::DINFO) << slice.as_string() << ", ";

    logger(LV::DINFO) << "\n";

    for (auto & train : simulator.get_trains())
    {
        for (auto & bunch : train.get_bunches())
        {
            auto const & ref = bunch.get_reference_particle();

            for (auto const & opn : operations)
            {
                logger(LV::INFO) << "    Independent_operator: operation type = " 
                    << opn->get_type() << std::endl;

                opn->apply(bunch, logger);

                simulator.diag_action_operation(*opn);
            }

            bunch.update_total_num();
        }
    }


#if 0
    if (do_update) 
    {
        if (verbosity > 3)
            logger << "Independent_operator: updating operations" << std::endl;

        update_operations(bunch.get_reference_particle());
        t = simple_timer_show(t, "independent_operator_apply-update_operations");
    }

    for (Independent_operations::iterator it = operations.begin();
            it != operations.end(); ++it) 
    {
        if (verbosity > 3)
            logger << "Independent_operator: operation type = "
                   << (*it)->get_type() << std::endl;

        (*it)->apply(bunch, verbosity, logger);

        std::string label("independent_operator_apply-" + (*it)->get_type()
                        + "_operation_apply");
        t = simple_timer_show(t, label.c_str());

        for (Diagnosticss::const_iterator itd =
                per_operation_diagnosticss.begin();
                itd != per_operation_diagnosticss.end(); ++it) 
        {
            (*itd)->update_and_write();
        }

        t = simple_timer_show(t, "diagnostics-operation");
    }

    bunch.update_total_num();
    t = simple_timer_show(t, "independent_operator_apply-bunch_update_total_num");
    t_total = simple_timer_show(t_total, "independent_operator_apply-total");
#endif
}

void Independent_operator::print_impl(Logger & logger) const
{
    logger(LoggerV::DEBUG)
        << "\tslices: " << "\n"
        << "\toperations: ";

    for(auto const & opn : operations) opn->print(logger);

    logger(LoggerV::DEBUG)
        << "\n";
}

#if 0
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
#endif

#if 0
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
#endif

#if 0
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
#endif

