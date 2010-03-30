#include "chef_lattice.h"
#include "chef_utils.h"
#include "utils/floating_point.h"
#include <beamline/beamline_elements.h>
#include <basic_toolkit/PhysicsConstants.h>
#include <physics_toolkit/DriftConverter.h>
#include <beamline/RefRegVisitor.h>

#include <stdexcept>

beamline
Chef_lattice::construct_raw_lattice(Lattice_element_to_chef_fn_map const& map)
{
    beamline raw_beamline;
    for (Lattice_elements::const_iterator latt_it =
            lattice_ptr->get_elements().begin(); latt_it
            != lattice_ptr->get_elements().end(); ++latt_it) {
        Lattice_element_to_chef_fn_map::const_iterator map_it = map.find(
                latt_it->get_type());
        if (map_it == map.end()) {
            throw(runtime_error("Chef_lattice: " + latt_it->get_type()
                    + " not handled"));
        } else {
            Chef_elements celms = map_it->second(*latt_it, brho);
            for (Chef_elements::const_iterator cel_it = celms.begin(); cel_it
                    != celms.end(); ++cel_it) {
                raw_beamline.append(*cel_it);
            }
            raw_beamline.append(lattice_element_marker);
        }
    }
    return raw_beamline;
}

void
Chef_lattice::polish_lattice(beamline const& raw_beamline)
{
    DriftConverter drift_converter;
    beamline_sptr = drift_converter.convert(raw_beamline);
    Particle testpart(reference_particle_to_chef_particle(
            lattice_ptr->get_reference_particle()));
    RefRegVisitor registrar(testpart);
    beamline_sptr->accept(registrar);
}

void
Chef_lattice::extract_element_map()
{
    Chef_elements chef_elements;
    Lattice_elements::iterator le_it = lattice_ptr->get_elements().begin();
    for (beamline::const_iterator b_it = beamline_sptr->begin(); b_it
            != beamline_sptr->end(); ++b_it) {
        if ((*b_it)->Name() == lattice_element_marker->Name()) {
            element_map[&(*le_it)] = chef_elements;
            chef_elements.clear();
            ++le_it;
        } else {
            chef_elements.push_back(*b_it);
            std::cout << "jfa: added chef_element " << (*b_it)->Name()
                    << " to lattice_element " << le_it->get_name() << std::endl;
        }
    }
}

void
Chef_lattice::construct(Lattice_element_to_chef_fn_map const& map)
{
    sliced_beamline_sptr = BmlPtr(new beamline("sliced"));
    if (!lattice_ptr->has_reference_particle()) {
        throw(std::runtime_error(
                "Chef_lattice: requires a reference particle in Lattice"));
    }
    brho = lattice_ptr->get_reference_particle().get_momentum()
            / PH_CNV_brho_to_p;

    polish_lattice(construct_raw_lattice(map));
    extract_element_map();
}

Chef_lattice::Chef_lattice(Lattice & lattice) :
    lattice_ptr(&lattice), beamline_sptr(), lattice_element_marker(new marker(
            "synergia_lattice_element_marker"))
{
    construct(get_standard_lattice_element_to_chef_fn_map());
}

Chef_lattice::Chef_lattice(Lattice & lattice,
        Lattice_element_to_chef_fn_map const& map) :
    lattice_ptr(&lattice), beamline_sptr(), lattice_element_marker(new marker(
            "synergia_lattice_element_marker"))
{
    construct(map);
}

Chef_elements &
Chef_lattice::get_chef_elements(Lattice_element const& lattice_element)
{
    return element_map[&lattice_element];
}

ElmPtr
slice_chef_element(ElmPtr & elm, double left, double right, double tolerance)
{
    std::cout << "jfa: slice_chef_element " << elm->Name() << " " << left
            << " " << right << std::endl;
    double length = elm->Length();
    ElmPtr retval, left_part, right_part;
    if (left == 0.0) {
        if (floating_point_equal(length, right, tolerance)) {
            retval = elm;
        } else {
            elm->Split(right / length, left_part, right_part);
            retval = left_part;
        }
    } else {
        elm->Split(left / length, left_part, right_part);
        if (floating_point_equal(length, right, tolerance)) {
            retval = right_part;
        } else {
            ElmPtr second_left_part, second_right_part;
            right_part->Split((right - left) / (length - left),
                    second_left_part, second_right_part);
            retval = second_left_part;
        }
    }

    return retval;
}

Chef_elements
Chef_lattice::get_chef_elements_from_slice(Lattice_element_slice & slice)
{
    Chef_elements all_elements = element_map[&(slice.get_lattice_element())];
    Chef_elements retval;
    if (slice.is_whole()) {
        std::cout << "jfa: slice is whole!\n";
        retval = all_elements;
    } else {
        const double tolerance = 1.0e-8;
        double left = slice.get_left();
        double right = slice.get_right();
        std::cout << "jfa: left,right = " << left << "," << right << std::endl;
        double s = 0.0;
        Chef_elements::iterator c_it = all_elements.begin();
        bool complete = false;
        double element_left, element_right;
        double total_done = 0.0;
        while (!complete) {
            std::cout << "jfa: working on " << (*c_it)->Name() << std::endl;
            double chef_element_length = (*c_it)->Length();
            if (!floating_point_leq(left, s + chef_element_length, tolerance)) {
                std::cout << "jfa !leq:" << s << " " << chef_element_length
                        << " " << left << " " << right << " " << std::endl;
                s += chef_element_length;
                ++c_it;
                if (c_it == all_elements.end()) {
                    throw(std::runtime_error(
                            "get_chef_elements_from_slice iterated beyond end of element list"));
                }
            } else {
                element_left = left - s;
                if (floating_point_leq(right, s + chef_element_length,
                        tolerance)) {
                    element_right = right - s;
                    retval.push_back(slice_chef_element(*c_it, element_left,
                            element_right, tolerance));
                    std::cout << "jfa: pushed back1 " << retval.back()->Name()
                            << std::endl;
                    total_done += element_right - element_left;
                } else {
                    element_right = chef_element_length;
                    retval.push_back(slice_chef_element(*c_it, element_left,
                            element_right, tolerance));
                    std::cout << "jfa: pushed back3 " << retval.back()->Name()
                            << std::endl;

                    s += chef_element_length;
                    ++c_it;
                }
                if (floating_point_equal(element_right, chef_element_length,
                        tolerance)) {
                    while ((++c_it != all_elements.end()) && ((*c_it)->Length()
                            == 0.0)) {
                        retval.push_back(*c_it);
                        std::cout << "jfa: pushed back2 "
                                << retval.back()->Name() << std::endl;

                    }
                }
                if (floating_point_equal(total_done, right - left, tolerance)) {
                    complete = true;
                }

            }
        }
    }

    std::cout << "jfa: end of get_chef_elements_from_slice\n";
    return retval;
}

void
Chef_lattice::construct_sliced_beamline(Lattice_element_slices & slices)
{
    std::cout << "jfa1 size =" << slices.size() << "\n";
    sliced_beamline_sptr->clear();
    for (Lattice_element_slices::iterator it = slices.begin(); it
            != slices.end(); ++it) {
        std::cout << "jfa2\n";
        Chef_elements chef_elements = get_chef_elements_from_slice(*(*it));
        for (Chef_elements::const_iterator c_it = chef_elements.begin(); c_it
                != chef_elements.end(); ++c_it) {
            sliced_beamline_sptr->append(*c_it);
        }
    }
}

BmlPtr
Chef_lattice::get_beamline_sptr()
{
    return beamline_sptr;
}

BmlPtr
Chef_lattice::get_sliced_beamline_sptr()
{
    return sliced_beamline_sptr;
}

Chef_lattice::~Chef_lattice()
{

}

Lattice_element_to_chef_fn_map
get_standard_lattice_element_to_chef_fn_map()
{
    Lattice_element_to_chef_fn_map map;
    map["marker"] = lattice_element_to_chef_marker;
    map["drift"] = lattice_element_to_chef_drift;
    map["sbend"] = lattice_element_to_chef_sbend;
    map["rbend"] = lattice_element_to_chef_rbend;
    map["quadrupole"] = lattice_element_to_chef_quadrupole;
    //    map["sextupole"] = lattice_element_to_chef_sextupole;
    //    map["octupole"] = lattice_element_to_chef_octupole;
    //    map["multipole"] = lattice_element_to_chef_multipole;
    //    map["solenoid"] = lattice_element_to_chef_solenoid;
    //    map["hkicker"] = lattice_element_to_chef_hkicker;
    //    map["vkicker"] = lattice_element_to_chef_vkicker;
    //    map["kicker"] = lattice_element_to_chef_kicker;
    //    map["rfcavity"] = lattice_element_to_chef_rfcavity;
    //    map["elseparator"] = lattice_element_to_chef_elseperator;
    //    map["hmonitor"] = lattice_element_to_chef_hmonitor;
    //    map["vmonitor"] = lattice_element_to_chef_vmonitor;
    //    map["monitor"] = lattice_element_to_chef_monitor;
    //    map["instrument"] = lattice_element_to_chef_instrument;
    //    map["ecollimator"] = lattice_element_to_chef_ecollimator;
    //    map["rcollimator"] = lattice_element_to_chef_rcollimator;

    return map;
}

Chef_elements
lattice_element_to_chef_marker(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(new marker(lattice_element.get_name().c_str()));
    retval.push_back(elm);
    return retval;
}

Chef_elements
lattice_element_to_chef_drift(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(new drift(lattice_element.get_name().c_str(),
            lattice_element.get_length()));
    retval.push_back(elm);
    return retval;
}

Chef_elements
lattice_element_to_chef_quadrupole(Lattice_element const& lattice_element,
        double brho)
{
    // tilt can have a string value of ""
    if (lattice_element.has_string_attribute("tilt")) {
        throw(runtime_error(
                "lattice_element_to_chef_quadrupole: tilt element not handled"));
    }
    if (lattice_element.has_double_attribute("tilt")) {
        if (lattice_element.get_double_attribute("tilt") != 0.0) {
            throw(runtime_error(
                    "lattice_element_to_chef_quadrupole: non-zero tilt element not handled"));
        }
    }
    Chef_elements retval;

    double length = lattice_element.get_length();
    bmlnElmnt* bmln_elmnt;
    if (length == 0.0) {
        bmln_elmnt = new thinQuad(lattice_element.get_name().c_str(), brho
                * lattice_element.get_double_attribute("k1"));
    } else {
        bmln_elmnt = new quadrupole(lattice_element.get_name().c_str(), length,
                brho * lattice_element.get_double_attribute("k1"));
    }
    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

Chef_elements
lattice_element_to_chef_sbend(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");

    if ((lattice_element.get_double_attribute("k1") != 0.0)
            || (lattice_element.get_double_attribute("k2") != 0.0)
            || (lattice_element.get_double_attribute("k3") != 0.0)
            || (lattice_element.get_double_attribute("tilt") != 0.0)
            || (lattice_element.get_double_attribute("h1") != 0.0)
            || (lattice_element.get_double_attribute("h2") != 0.0)
            || (lattice_element.get_double_attribute("hgap") != 0.0)
            || (lattice_element.get_double_attribute("fint") != 0.0)) {
        throw(runtime_error(
                "lattice_element_to_chef_sbend: non-zero element(s) of something not handled"));
    }

    ElmPtr elm(new sbend(lattice_element.get_name().c_str(), length, brho
            * angle / length, angle, e1, e2));
    retval.push_back(elm);
    return retval;
}

Chef_elements
lattice_element_to_chef_rbend(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");

    if ((lattice_element.get_double_attribute("k1") != 0.0)
            || (lattice_element.get_double_attribute("k2") != 0.0)
            || (lattice_element.get_double_attribute("k3") != 0.0)
            || (lattice_element.get_double_attribute("tilt") != 0.0)
            || (lattice_element.get_double_attribute("h1") != 0.0)
            || (lattice_element.get_double_attribute("h2") != 0.0)
            || (lattice_element.get_double_attribute("hgap") != 0.0)
            || (lattice_element.get_double_attribute("fint") != 0.0)) {
        throw(runtime_error(
                "lattice_element_to_chef_sbend: non-zero element(s) of something not handled"));
    }

    ElmPtr elm(new rbend(lattice_element.get_name().c_str(), length, brho
            * (2.0 * sin(0.5 * angle)) / length, angle, e1, e2));
    retval.push_back(elm);
    return retval;
}

