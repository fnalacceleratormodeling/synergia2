#include "chef_lattice.h"
#include "chef_utils.h"
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/math_constants.h"
#include <beamline/beamline_elements.h>
#include <basic_toolkit/PhysicsConstants.h>
#include <physics_toolkit/DriftConverter.h>
#include <beamline/RefRegVisitor.h>

#include <stdexcept>

void
Chef_lattice::construct_beamline()
{
    BmlPtr unpolished_beamline_sptr(new beamline());
    for (Lattice_elements::const_iterator latt_it =
            lattice_sptr->get_elements().begin(); latt_it
            != lattice_sptr->get_elements().end(); ++latt_it) {
        std::string type((*latt_it)->get_type());
        if (!element_adaptor_map_sptr->has_adaptor(type)) {
            throw(runtime_error("Chef_lattice: " + type + " not handled"));
        } else {
            Chef_elements
                    celms =
                            element_adaptor_map_sptr->get_adaptor(type)->get_chef_elements(
                                    *(*latt_it), brho);
            for (Chef_elements::const_iterator cel_it = celms.begin(); cel_it
                    != celms.end(); ++cel_it) {
                unpolished_beamline_sptr->append(*cel_it);
            }
            unpolished_beamline_sptr->append(lattice_element_marker);
        }
    }
    beamline_sptr = polish_beamline(unpolished_beamline_sptr);
}

void
Chef_lattice::register_beamline(beamline & the_beamline)
{
    Particle testpart(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    RefRegVisitor registrar(testpart);
    the_beamline.accept(registrar);
}

BmlPtr
Chef_lattice::polish_beamline(BmlPtr beamline_sptr)
{
    DriftConverter drift_converter;
    BmlPtr converted_beamline_sptr;
    if (beamline_sptr->countHowManyDeeply() < 3) {
        converted_beamline_sptr = beamline_sptr;
    } else {
        converted_beamline_sptr = drift_converter.convert(*beamline_sptr);
        register_beamline(*converted_beamline_sptr);
    }
    return converted_beamline_sptr;
}

void
Chef_lattice::extract_element_map()
{
    Chef_elements chef_elements;
    Lattice_elements::iterator le_it = lattice_sptr->get_elements().begin();
    for (beamline::const_iterator b_it = beamline_sptr->begin(); b_it
            != beamline_sptr->end(); ++b_it) {
        if ((*b_it)->Name() == lattice_element_marker->Name()) {
            element_map[le_it->get()] = chef_elements;
            chef_elements.clear();
            ++le_it;
        } else {
            chef_elements.push_back(*b_it);
        }
    }
}

void
Chef_lattice::extract_element_slice_map(Lattice_element_slices const& slices)
{
    Chef_elements chef_elements;
    Lattice_element_slices::const_iterator s_it = slices.begin();
    for (beamline::const_iterator b_it = sliced_beamline_sptr->begin(); b_it
            != sliced_beamline_sptr->end(); ++b_it) {
        if ((*b_it)->Name() == lattice_element_marker->Name()) {
            element_slice_map[s_it->get()] = chef_elements;
            chef_elements.clear();
            ++s_it;
        } else {
            chef_elements.push_back(*b_it);
        }
    }
}

void
Chef_lattice::construct()
{
    lattice_sptr->complete_attributes(*element_adaptor_map_sptr);
    sliced_beamline_sptr = BmlPtr(new beamline("sliced"));
    have_sliced_beamline_ = false;
    if (!lattice_sptr->has_reference_particle()) {
        throw(std::runtime_error(
                "Chef_lattice: requires a reference particle in Lattice"));
    }
    brho = lattice_sptr->get_reference_particle().get_momentum()
            / PH_CNV_brho_to_p;

    construct_beamline();
    extract_element_map();
}

Chef_lattice::Chef_lattice(Lattice_sptr lattice_sptr) :
            lattice_sptr(lattice_sptr),
            beamline_sptr(),
            lattice_element_marker(
                    new marker("synergia_lattice_element_marker")),
            element_adaptor_map_sptr(new Element_adaptor_map)
{
    construct();
}

Chef_lattice::Chef_lattice(Lattice_sptr lattice_sptr,
        Element_adaptor_map_sptr element_adaptor_map_sptr) :
            lattice_sptr(lattice_sptr),
            beamline_sptr(),
            lattice_element_marker(
                    new marker("synergia_lattice_element_marker")),
            element_adaptor_map_sptr(element_adaptor_map_sptr)
{
    construct();
}

Element_adaptor_map_sptr
Chef_lattice::get_element_adaptor_map_sptr()
{
    return element_adaptor_map_sptr;
}

Chef_elements &
Chef_lattice::get_chef_elements(Lattice_element const& lattice_element)
{
    return element_map[&lattice_element];
}

Chef_elements &
Chef_lattice::get_chef_elements(
        Lattice_element_slice const& lattice_element_slice)
{
    if (!have_sliced_beamline_) {
        throw std::runtime_error(
                "get_chef_elements(Lattice_element_slice const&) called before construct_sliced_beamline\n");
    }
    if (element_slice_map.count(&lattice_element_slice) == 0) {
        throw std::runtime_error(
                "get_chef_elements(Lattice_element_slice const&): slice not found\n");
    }
    return element_slice_map[&lattice_element_slice];
}

ElmPtr
slice_chef_element(ElmPtr & elm, double left, double right, double tolerance)
{
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
Chef_lattice::get_chef_elements_from_slice(Lattice_element_slice const& slice)
{
    Chef_elements all_elements = element_adaptor_map_sptr->get_adaptor(
            slice.get_lattice_element().get_type())->get_chef_elements(
            slice.get_lattice_element(), brho);
    Chef_elements retval;
    if (slice.is_whole()) {
        retval = all_elements;
    } else {
        const double tolerance = 1.0e-8;
        double left = slice.get_left();
        double right = slice.get_right();
        double s = 0.0;
        Chef_elements::iterator c_it = all_elements.begin();
        bool complete = false;
        double element_left, element_right;
        double total_done = 0.0;
        while (!complete) {
            element_left = left - s + total_done;
            double chef_element_length = (*c_it)->Length();
            if (floating_point_leq(chef_element_length, element_left, tolerance)
                    && (total_done == 0.0)) {
                s += chef_element_length;
                ++c_it;
                if (c_it == all_elements.end()) {
                    throw(std::runtime_error(
                            "get_chef_elements_from_slice iterated beyond end of element list"));
                }
            } else {
                if (chef_element_length == 0.0) {
                    retval.push_back(*c_it);
                    ++c_it;
                } else {
                    if (floating_point_leq(right, s + chef_element_length,
                            tolerance)) {
                        // take part of element
                        element_right = right - s;
                        retval.push_back(
                                slice_chef_element(*c_it, element_left,
                                        element_right, tolerance));
                        s += element_right - element_left;
                        total_done += element_right - element_left;
                    } else {
                        // take rest of element
                        element_right = chef_element_length;
                        retval.push_back(
                                slice_chef_element(*c_it, element_left,
                                        element_right, tolerance));
                        s += chef_element_length;
                        total_done += element_right - element_left;
                        ++c_it;
                    }
                }
                if (floating_point_equal(total_done, right - left, tolerance)) {
                    complete = true;
                    if (floating_point_equal(element_right,
                            chef_element_length, tolerance)) {
                        while ((++c_it != all_elements.end())
                                && ((*c_it)->Length() == 0.0)) {
                            retval.push_back(*c_it);
                        }
                    }
                }
            }
        }
    }

    return retval;
}

bool
Chef_lattice::have_sliced_beamline() const
{
    return have_sliced_beamline_;
}

void
Chef_lattice::construct_sliced_beamline(Lattice_element_slices const& slices)
{
    BmlPtr unpolished_beamline_sptr(new beamline());
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        Chef_elements chef_elements = get_chef_elements_from_slice(*(*it));
        for (Chef_elements::const_iterator c_it = chef_elements.begin(); c_it
                != chef_elements.end(); ++c_it) {
            unpolished_beamline_sptr->append(*c_it);
        }
        unpolished_beamline_sptr->append(lattice_element_marker);
    }
    sliced_beamline_sptr = polish_beamline(unpolished_beamline_sptr);
    extract_element_slice_map(slices);
    have_sliced_beamline_ = true;
}

BmlPtr
Chef_lattice::get_beamline_sptr()
{
    return beamline_sptr;
}

BmlPtr
Chef_lattice::get_sliced_beamline_sptr()
{
    if (!have_sliced_beamline_) {
        throw std::runtime_error(
                "get_sliced_beamline_sptr() called before construct_sliced_beamline\n");
    }
    return sliced_beamline_sptr;
}

void
print_chef_beamline(BmlPtr const& chef_beamline,
        Reference_particle const& reference_particle)
{
    Particle testpart(reference_particle_to_chef_particle(reference_particle));
    for (beamline::deep_iterator it = chef_beamline->deep_begin(); it
            != chef_beamline->deep_end(); ++it) {
        std::cout << (*it)->Name() << ": " << (*it)->Type() << ", length = "
                << (*it)->OrbitLength(testpart) << ", strength = "
                << (*it)->Strength() << std::endl;
        (*it)->propagate(testpart);
    }
}

void
Chef_lattice::print_beamline() const
{
    print_chef_beamline(beamline_sptr, lattice_sptr->get_reference_particle());
}

void
Chef_lattice::print_sliced_beamline() const
{
    print_chef_beamline(sliced_beamline_sptr,
            lattice_sptr->get_reference_particle());
}

Chef_lattice::~Chef_lattice()
{

}
