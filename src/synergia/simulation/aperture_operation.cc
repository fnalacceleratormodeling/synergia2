#include "aperture_operation.h"
#include "synergia/foundation/math_constants.h"
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <cmath>

Aperture_operation::Aperture_operation(Lattice_element_slice_sptr slice) :
    Independent_operation("aperture"), slice(slice), deposited_charge(0.0)
{
}

void
Aperture_operation::reset_deposited_charge()
{
    deposited_charge = 0.0;
}

void
Aperture_operation::increase_deposited_charge(double charge)
{
    deposited_charge += charge;
}

double
Aperture_operation::get_deposited_charge() const
{
    return deposited_charge;
}

Aperture_operation::~Aperture_operation()
{
}

const char Finite_aperture_operation::aperture_type[] = "finite";

Finite_aperture_operation::Finite_aperture_operation(
        Lattice_element_slice_sptr slice) :
    Aperture_operation(slice)
{
}

const char *
Finite_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Finite_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    return (aperture_type == aperture_operation.get_aperture_type());
}

void
Finite_aperture_operation::apply(Bunch & bunch)
{
    apply_impl(*this, bunch);
}

Finite_aperture_operation::~Finite_aperture_operation()
{
}

const double Circular_aperture_operation::default_radius(1000.0);
const char Circular_aperture_operation::aperture_type[] = "circular";
const char Circular_aperture_operation::attribute_name[] = "circular";

Circular_aperture_operation::Circular_aperture_operation(
        Lattice_element_slice_sptr slice) :
    Aperture_operation(slice)
{

    if (slice->get_lattice_element().has_double_attribute(
            "circular_aperture_radius")) {
        radius = slice->get_lattice_element().get_double_attribute(
                "circular_aperture_radius");
    } else {
        radius = default_radius;
    }
    radius2 = radius * radius;
}

const char *
Circular_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Circular_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Circular_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Circular_aperture_operation::operator==(
        Circular_aperture_operation const& circular_aperture_operation) const
{
    return (radius == circular_aperture_operation.radius);
}

void
Circular_aperture_operation::apply(Bunch & bunch)
{
    apply_impl(*this, bunch);
}

Circular_aperture_operation::~Circular_aperture_operation()
{
}

const char Elliptical_aperture_operation::aperture_type[] = "elliptical";
const char Elliptical_aperture_operation::attribute_name[] = "elliptical";

Elliptical_aperture_operation::Elliptical_aperture_operation(
        Lattice_element_slice_sptr slice) :
    Aperture_operation(slice)
{
    if (slice->get_lattice_element().has_double_attribute(
            "elliptical_aperture_horizontal_radius")) {
        horizontal_radius = slice->get_lattice_element().get_double_attribute(
                "elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_horizontal_radius attribute");
    }
    if (slice->get_lattice_element().has_double_attribute(
            "elliptical_aperture_vertical_radius")) {
        vertical_radius = slice->get_lattice_element().get_double_attribute(
                "elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_vertical_radius attribute");
    }
    h2 = horizontal_radius * horizontal_radius;
    v2 = vertical_radius * vertical_radius;
}

const char *
Elliptical_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Elliptical_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Elliptical_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Elliptical_aperture_operation::operator==(
        Elliptical_aperture_operation const& elliptical_aperture_operation) const
{
    return ((horizontal_radius
            == elliptical_aperture_operation.horizontal_radius)
            && (vertical_radius
                    == elliptical_aperture_operation.vertical_radius));
}

void
Elliptical_aperture_operation::apply(Bunch & bunch)
{
    apply_impl(*this, bunch);
}

Elliptical_aperture_operation::~Elliptical_aperture_operation()
{
}

const char Rectangular_aperture_operation::aperture_type[] = "rectangular";
const char Rectangular_aperture_operation::attribute_name[] = "rectangular";

Rectangular_aperture_operation::Rectangular_aperture_operation(
        Lattice_element_slice_sptr slice) :
    Aperture_operation(slice)
{
    if (slice->get_lattice_element().has_double_attribute(
            "rectangular_aperture_width")) {
        width = slice->get_lattice_element().get_double_attribute(
                "rectangular_aperture_width");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_width attribute");
    }
    if (slice->get_lattice_element().has_double_attribute(
            "rectangular_aperture_height")) {
        height = slice->get_lattice_element().get_double_attribute(
                "rectangular_aperture_height");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_height attribute");
    }
}

const char *
Rectangular_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Rectangular_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Rectangular_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Rectangular_aperture_operation::operator==(
        Rectangular_aperture_operation const& rectangular_aperture_operation) const
{
    return ((width == rectangular_aperture_operation.width) && (height
            == rectangular_aperture_operation.height));
}

void
Rectangular_aperture_operation::apply(Bunch & bunch)
{
    apply_impl(*this, bunch);
}

Rectangular_aperture_operation::~Rectangular_aperture_operation()
{
}

const char Polygon_aperture_operation::aperture_type[] = "polygon";
const char Polygon_aperture_operation::attribute_name[] = "polygon";

Polygon_aperture_operation::Polygon_aperture_operation(
        Lattice_element_slice_sptr slice) :
    Aperture_operation(slice)
{
    if (slice->get_lattice_element().has_double_attribute(
            "the_number_of_vertices")) {
        vertices_num = int(
                slice->get_lattice_element().get_double_attribute(
                        "the_number_of_vertices"));
        if (vertices_num < 3) throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires at least 3 vertices");
    } else {
        throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires the_number_of vertices attribute");
    }
    for (int index = 0; index < vertices_num; ++index) {
        std::string ss = boost::lexical_cast<std::string >(index + 1);
        std::string x = "pax" + ss;
        std::string y = "pay" + ss;
        if ((slice->get_lattice_element().has_double_attribute(x))
                && (slice->get_lattice_element().has_double_attribute(y))) {
            vertices.push_back(
                    std::complex<double >(
                            slice->get_lattice_element().get_double_attribute(x),
                            slice->get_lattice_element().get_double_attribute(y)));
            //std::cout << index << "  " << slice->get_lattice_element().get_double_attribute(x)
            //        << "  " << slice->get_lattice_element().get_double_attribute(y) << std::endl;
        } else {
            throw std::runtime_error(
                    "Polygon_aperture_operation: polygon_aperture requires x and y coordinate attributes for each vertex");
        }
    }
}

const char *
Polygon_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Polygon_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Polygon_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Polygon_aperture_operation::operator==(
        Polygon_aperture_operation const& polygon_aperture_operation) const
{
    return (vertices == polygon_aperture_operation.vertices);
}

void
Polygon_aperture_operation::apply(Bunch & bunch)
{
    apply_impl(*this, bunch);
}

Polygon_aperture_operation::~Polygon_aperture_operation()
{
}

const char Wire_elliptical_aperture_operation::aperture_type[] =
        "wire_elliptical";
const char Wire_elliptical_aperture_operation::attribute_name[] =
        "wire_elliptical";

Wire_elliptical_aperture_operation::Wire_elliptical_aperture_operation(
        Lattice_element_slice_sptr slice) :
    Aperture_operation(slice)
{
    if (slice->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_horizontal_radius")) {
        horizontal_radius = slice->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_horizontal_radius attribute");
    }
    if (slice->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_vertical_radius")) {
        vertical_radius = slice->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_vertical_radius attribute");
    }
    if (slice->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_wire_x")) {
        wire_x = slice->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_wire_x");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_x attribute");
    }
    if (slice->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_wire_width")) {
        wire_width = slice->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_wire_width");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_width attribute");
    }
    if (slice->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_gap")) {
        gap = slice->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_gap");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_gap attribute");
    }
    h2 = horizontal_radius * horizontal_radius;
    v2 = vertical_radius * vertical_radius;
}

const char *
Wire_elliptical_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Wire_elliptical_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Wire_elliptical_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Wire_elliptical_aperture_operation::operator==(
        Wire_elliptical_aperture_operation const& wire_elliptical_aperture_operation) const
{
    return ((horizontal_radius
            == wire_elliptical_aperture_operation.horizontal_radius)
            && (vertical_radius
                    == wire_elliptical_aperture_operation.vertical_radius)
            && (wire_x == wire_elliptical_aperture_operation.wire_x)
            && (wire_width == wire_elliptical_aperture_operation.wire_width)
            && (gap == wire_elliptical_aperture_operation.gap));
}

void
Wire_elliptical_aperture_operation::apply(Bunch & bunch)
{
    apply_impl(*this, bunch);
}

Wire_elliptical_aperture_operation::~Wire_elliptical_aperture_operation()
{
}
