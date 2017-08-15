#include "aperture_operation.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/simple_timer.h"
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <cmath>

const char Aperture_operation::charge_attribute[] = "deposited_charge";

Aperture_operation::Aperture_operation(Lattice_element_slice_sptr slice_sptr) :
    Independent_operation("aperture"), slice_sptr(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute("hoffset")) {
        x_offset = slice_sptr->get_lattice_element().get_double_attribute("hoffset");
    } else {
        x_offset = 0.0;
    }
    if (slice_sptr->get_lattice_element().has_double_attribute("voffset")) {
        y_offset = slice_sptr->get_lattice_element().get_double_attribute("voffset");
    } else {
        y_offset = 0.0;
    }
}

Aperture_operation::Aperture_operation() :
    Independent_operation("aperture")
{
}

Lattice_element_slice_sptr 
Aperture_operation::get_slice_sptr() const
{
  return slice_sptr;
} 

void
Aperture_operation::deposit_charge(double charge)
{
    double deposited_charge(0.0);
    if (slice_sptr->get_lattice_element().has_double_attribute(charge_attribute)) {
        deposited_charge
                = slice_sptr->get_lattice_element().get_double_attribute(
                        charge_attribute);
    }
    deposited_charge += charge;
    slice_sptr->get_lattice_element().set_double_attribute(charge_attribute,
            deposited_charge, false);
}

template<class Archive>
    void
    Aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Independent_operation);
        ar & BOOST_SERIALIZATION_NVP(slice_sptr);
        ar & BOOST_SERIALIZATION_NVP(x_offset);
        ar & BOOST_SERIALIZATION_NVP(y_offset);
    }

template
void
Aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Aperture_operation::~Aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Aperture_operation);

const char Finite_aperture_operation::aperture_type[] = "finite";

Finite_aperture_operation::Finite_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
}

Finite_aperture_operation::Finite_aperture_operation()
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
Finite_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "finite_aperture-apply");
}

template<class Archive>
    void
    Finite_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
    }

template
void
Finite_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Finite_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Finite_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Finite_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Finite_aperture_operation::~Finite_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Finite_aperture_operation);

const double Circular_aperture_operation::default_radius(1000.0);
const char Circular_aperture_operation::aperture_type[] = "circular";
const char Circular_aperture_operation::attribute_name[] = "circular";

Circular_aperture_operation::Circular_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{

    if (slice_sptr->get_lattice_element().has_double_attribute(
            "circular_aperture_radius")) {
        radius = slice_sptr->get_lattice_element().get_double_attribute(
                "circular_aperture_radius");
    } else {
        radius = default_radius;
    }
    radius2 = radius * radius;
}

Circular_aperture_operation::Circular_aperture_operation()
{
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
Circular_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "circular_aperture-apply");
}

template<class Archive>
    void
    Circular_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(radius);
        ar & BOOST_SERIALIZATION_NVP(radius2);
    }

template
void
Circular_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Circular_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Circular_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Circular_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Circular_aperture_operation::~Circular_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Circular_aperture_operation);

const char Elliptical_aperture_operation::aperture_type[] = "elliptical";
const char Elliptical_aperture_operation::attribute_name[] = "elliptical";

Elliptical_aperture_operation::Elliptical_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "elliptical_aperture_horizontal_radius")) {
        horizontal_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_horizontal_radius attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "elliptical_aperture_vertical_radius")) {
        vertical_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_vertical_radius attribute");
    }
    h2 = horizontal_radius * horizontal_radius;
    v2 = vertical_radius * vertical_radius;
}

Elliptical_aperture_operation::Elliptical_aperture_operation()
{
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
Elliptical_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "elliptical_aperture-apply");
}

template<class Archive>
    void
    Elliptical_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(horizontal_radius);
        ar & BOOST_SERIALIZATION_NVP(vertical_radius);
        ar & BOOST_SERIALIZATION_NVP(h2);
        ar & BOOST_SERIALIZATION_NVP(v2);
    }

template
void
Elliptical_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Elliptical_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Elliptical_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Elliptical_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Elliptical_aperture_operation::~Elliptical_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Elliptical_aperture_operation);

const char Rectangular_aperture_operation::aperture_type[] = "rectangular";
const char Rectangular_aperture_operation::attribute_name[] = "rectangular";

Rectangular_aperture_operation::Rectangular_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_width")) {
        width = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_width");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_width attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_height")) {
        height = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_height");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_height attribute");
    }
}

Rectangular_aperture_operation::Rectangular_aperture_operation()
{
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
Rectangular_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "rectangular_aperture-apply");
}

template<class Archive>
    void
    Rectangular_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(width);
        ar & BOOST_SERIALIZATION_NVP(height);
    }

template
void
Rectangular_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rectangular_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rectangular_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rectangular_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rectangular_aperture_operation::~Rectangular_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rectangular_aperture_operation);

const char Rectangular_with_ears_aperture_operation::aperture_type[] = "rectangular_with_ears";
const char Rectangular_with_ears_aperture_operation::attribute_name[] = "rectangular_with_ears";

Rectangular_with_ears_aperture_operation::Rectangular_with_ears_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_width")) {
        width = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_width");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_width attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_height")) {
        height = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_height");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_height attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_ear_offset")) {
        ear_offset = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_ear_offset");
    } else {
        throw std::runtime_error(
                "Rectangular_with_ears_aperture_operation: rectangular_with_ears_aperture requires an rectangular_aperture_ear_offset attribute");
    }
    radius = 0.5*height - ear_offset;
}

Rectangular_with_ears_aperture_operation::Rectangular_with_ears_aperture_operation()
{
}

const char *
Rectangular_with_ears_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Rectangular_with_ears_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Rectangular_with_ears_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Rectangular_with_ears_aperture_operation::operator==(
        Rectangular_with_ears_aperture_operation const& rectangular_with_ears_aperture_operation) const
{
    return ((width == rectangular_with_ears_aperture_operation.width) &&
            (height == rectangular_with_ears_aperture_operation.height) &&
            (ear_offset == rectangular_with_ears_aperture_operation.ear_offset));
}

void
Rectangular_with_ears_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "rectangular_aperture-apply");
}

template<class Archive>
    void
    Rectangular_with_ears_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(width);
        ar & BOOST_SERIALIZATION_NVP(height);
        ar & BOOST_SERIALIZATION_NVP(ear_offset);
        ar & BOOST_SERIALIZATION_NVP(radius);
    }

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rectangular_with_ears_aperture_operation::~Rectangular_with_ears_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rectangular_with_ears_aperture_operation);

const char Polygon_aperture_operation::aperture_type[] = "polygon";
const char Polygon_aperture_operation::attribute_name[] = "polygon";

Polygon_aperture_operation::Polygon_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "the_number_of_vertices")) {
        num_vertices = int(
                slice_sptr->get_lattice_element().get_double_attribute(
                        "the_number_of_vertices"));
        if (num_vertices < 3) throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires at least 3 vertices");
    } else {
        throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires the_number_of vertices attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "min_radius2")) {
        min_radius2 = slice_sptr->get_lattice_element().get_double_attribute(
                        "min_radius2");
    } else {
        min_radius2 = 0.0;
    }
    for (int index = 0; index < num_vertices; ++index) {
        std::string ss = boost::lexical_cast<std::string >(index + 1);
        std::string x = "pax" + ss;
        std::string y = "pay" + ss;
        if ((slice_sptr->get_lattice_element().has_double_attribute(x))
                && (slice_sptr->get_lattice_element().has_double_attribute(y))) {
            vertices.push_back(
                    std::complex<double >(
                            slice_sptr->get_lattice_element().get_double_attribute(
                                    x),
                            slice_sptr->get_lattice_element().get_double_attribute(
                                    y)));
        } else {
            throw std::runtime_error(
                    "Polygon_aperture_operation: polygon_aperture requires x and y coordinate attributes for each vertex");
        }
    }
}

Polygon_aperture_operation::Polygon_aperture_operation()
{
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
    return ((vertices == polygon_aperture_operation.vertices)
            && (min_radius2 == polygon_aperture_operation.min_radius2));
}

void
Polygon_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "polygon_aperture-apply");
}

template<class Archive>
    void
    Polygon_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(num_vertices);
        ar & BOOST_SERIALIZATION_NVP(vertices);
        ar & BOOST_SERIALIZATION_NVP(min_radius2);
    }

template
void
Polygon_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Polygon_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Polygon_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Polygon_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Polygon_aperture_operation::~Polygon_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Polygon_aperture_operation);

const char Wire_elliptical_aperture_operation::aperture_type[] =
        "wire_elliptical";
const char Wire_elliptical_aperture_operation::attribute_name[] =
        "wire_elliptical";

Wire_elliptical_aperture_operation::Wire_elliptical_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_horizontal_radius")) {
        horizontal_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "wire_elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_horizontal_radius attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_vertical_radius")) {
        vertical_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "wire_elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_vertical_radius attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_wire_x")) {
        wire_x = slice_sptr->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_wire_x");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_x attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_wire_width")) {
        wire_width = slice_sptr->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_wire_width");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_width attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_gap")) {
        gap = slice_sptr->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_gap");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_gap attribute");
    }
    h2 = horizontal_radius * horizontal_radius;
    v2 = vertical_radius * vertical_radius;
}

Wire_elliptical_aperture_operation::Wire_elliptical_aperture_operation()
{
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
Wire_elliptical_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "wire_elliptical_aperture-apply");
}

template<class Archive>
    void
    Wire_elliptical_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(horizontal_radius);
        ar & BOOST_SERIALIZATION_NVP(vertical_radius);
        ar & BOOST_SERIALIZATION_NVP(h2);
        ar & BOOST_SERIALIZATION_NVP(v2);
        ar & BOOST_SERIALIZATION_NVP(wire_x);
        ar & BOOST_SERIALIZATION_NVP(wire_width);
        ar & BOOST_SERIALIZATION_NVP(gap);
    }

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Wire_elliptical_aperture_operation::~Wire_elliptical_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Wire_elliptical_aperture_operation);

const char Lambertson_aperture_operation::aperture_type[] = "lambertson";
const char Lambertson_aperture_operation::attribute_name[] = "lambertson";

Lambertson_aperture_operation::Lambertson_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "lambertson_aperture_radius")) {
        radius = slice_sptr->get_lattice_element().get_double_attribute(
                "lambertson_aperture_radius");
    } else {
        throw std::runtime_error(
                "lambertson_aperture_operation: lambertson_aperture requires an lambertson_aperture_radius attribute");
    }
}

Lambertson_aperture_operation::Lambertson_aperture_operation()
{
}

const char *
Lambertson_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Lambertson_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Lambertson_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Lambertson_aperture_operation::operator==(
        Lambertson_aperture_operation const& lambertson_aperture_operation) const
{
    return (radius == lambertson_aperture_operation.radius);
}

void
Lambertson_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    dump_particles(*this, bunch, verbosity, logger);
    //apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "lambertson_aperture-apply");
}

template<class Archive>
    void
    Lambertson_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(radius);
    }

template
void
Lambertson_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lambertson_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);
    
template
void
Lambertson_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lambertson_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Lambertson_aperture_operation::~Lambertson_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Lambertson_aperture_operation)
