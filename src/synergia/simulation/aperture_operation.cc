#include "aperture_operation.h"
#include "synergia/foundation/math_constants.h"
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <cmath>

Aperture_operation::Aperture_operation(Lattice_element const& element) :
    Independent_operation("aperture")
{
}

Aperture_operation::~Aperture_operation()
{
}

const double Circular_aperture_operation::default_radius(1000.0);
const char Circular_aperture_operation::type_name[] = "circular_aperture";
const char Circular_aperture_operation::attribute_name[] = "circular";

Circular_aperture_operation::Circular_aperture_operation(
        Lattice_element const& element) :
    Aperture_operation(element)
{
    if (element.has_double_attribute("circular_aperture_radius")) {
        radius = element.get_double_attribute("circular_aperture_radius");
    } else {
        radius = default_radius;
    }
}

const char *
Circular_aperture_operation::get_type_name() const
{
    return type_name;
}

bool
Circular_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (type_name == aperture_operation.get_type_name()) {
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
    double radius2 = radius * radius;
    MArray2d_ref particles(bunch.get_local_particles());
    int kept = 0;
    int discarded = 0;
    int local_num = bunch.get_local_num();
    for (int part = 0; part < local_num; ++part) {
        bool try_discard = true;
        while (try_discard) {
            double r2 = particles[part][Bunch::x] * particles[part][Bunch::x]
                    + particles[part][Bunch::y] * particles[part][Bunch::y];
            if (r2 > radius2) {
                ++discarded;
                --local_num;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else {
                ++kept;
                try_discard = false;
            }
        }
    }
    //std::cout << "kept = " << kept << ", discarded = " << discarded
    //        << std::endl;
    bunch.set_local_num(local_num);
    //std::cout<<"circular aperture applied"<<std::endl;
}

Circular_aperture_operation::~Circular_aperture_operation()
{
}

const char Elliptical_aperture_operation::type_name[] = "elliptical_aperture";
const char Elliptical_aperture_operation::attribute_name[] = "elliptical";

Elliptical_aperture_operation::Elliptical_aperture_operation(
        Lattice_element const& element) :
    Aperture_operation(element)
{
    if (element.has_double_attribute("elliptical_aperture_horizontal_radius")) {
        horizontal_radius = element.get_double_attribute(
                "elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_horizontal_radius attribute");
    }
    if (element.has_double_attribute("elliptical_aperture_vertical_radius")) {
        vertical_radius = element.get_double_attribute(
                "elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_vertical_radius attribute");
    }
}

const char *
Elliptical_aperture_operation::get_type_name() const
{
    return type_name;
}

bool
Elliptical_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (type_name == aperture_operation.get_type_name()) {
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
    double h2 = horizontal_radius * horizontal_radius;
    double v2 = vertical_radius * vertical_radius;
    MArray2d_ref particles(bunch.get_local_particles());
    int kept = 0;
    int discarded = 0;
    int local_num = bunch.get_local_num();
    for (int part = 0; part < local_num; ++part) {
        bool try_discard = true;
        while (try_discard) {
            double scaled_r2 = particles[part][Bunch::x]
                    * particles[part][Bunch::x] / h2
                    + particles[part][Bunch::y] * particles[part][Bunch::y]
                            / v2;
            if (scaled_r2 > 1.0) {
                ++discarded;
                --local_num;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else {
                ++kept;
                try_discard = false;
            }
        }
    }
    //std::cout << "kept = " << kept << ", discarded = " << discarded
    //        << std::endl;
    bunch.set_local_num(local_num);
}

Elliptical_aperture_operation::~Elliptical_aperture_operation()
{
}

const char Rectangular_aperture_operation::type_name[] = "rectangular_aperture";
const char Rectangular_aperture_operation::attribute_name[] = "rectangular";

Rectangular_aperture_operation::Rectangular_aperture_operation(
        Lattice_element const& element) :
    Aperture_operation(element)
{
    if (element.has_double_attribute("rectangular_aperture_width")) {
        width = element.get_double_attribute("rectangular_aperture_width");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_width attribute");
    }
    if (element.has_double_attribute("rectangular_aperture_height")) {
        height = element.get_double_attribute("rectangular_aperture_height");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_height attribute");
    }
}

const char *
Rectangular_aperture_operation::get_type_name() const
{
    return type_name;
}

bool
Rectangular_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (type_name == aperture_operation.get_type_name()) {
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
    MArray2d_ref particles(bunch.get_local_particles());
    int kept = 0;
    int discarded = 0;
    int local_num = bunch.get_local_num();
    for (int part = 0; part < local_num; ++part) {
        bool try_discard = true;
        while (try_discard) {
            if ((std::abs(particles[part][Bunch::x]) > 0.5 * width)
                    || (std::abs(particles[part][Bunch::y]) > 0.5 * height)) {
                ++discarded;
                --local_num;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else {
                ++kept;
                try_discard = false;
            }
        }
    }
    //std::cout << "kept = " << kept << ", discarded = " << discarded
    //        << std::endl;
    bunch.set_local_num(local_num);
    //std::cout<<"rectangular aperture applied"<<std::endl;
}

Rectangular_aperture_operation::~Rectangular_aperture_operation()
{
}

const char Polygon_aperture_operation::type_name[] = "polygon_aperture";
const char Polygon_aperture_operation::attribute_name[] = "polygon";

Polygon_aperture_operation::Polygon_aperture_operation(
        Lattice_element const& element) :
    Aperture_operation(element)
{
    if (element.has_double_attribute("the_number_of_vertices")) {
        vertices_num = int(element.get_double_attribute("the_number_of_vertices"));
        if (vertices_num < 3)
                throw std::runtime_error(
                        "Polygon_aperture_operation: polygon_aperture requires at least 3 vertices");
    } else {
        throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires the_number_of vertices attribute");
    }
    for (int index = 0; index < vertices_num; ++index) {
        std::string ss = boost::lexical_cast<std::string>(index + 1);
        std::string x = "pax" + ss;
        std::string y = "pay" + ss;
        if ((element.has_double_attribute(x)) && (element.has_double_attribute(y))) {
            vertices.push_back(std::complex<double > (element.get_double_attribute(x), element.get_double_attribute(y)));
            //std::cout << index << "  " << element.get_double_attribute(x)
            //        << "  " << element.get_double_attribute(y) << std::endl;
        } else {
            throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires x and y coordinate attributes for each vertex");
        }
    }
}

const char *
Polygon_aperture_operation::get_type_name() const
{
    return type_name;
}

bool
Polygon_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (type_name == aperture_operation.get_type_name()) {
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
    MArray2d_ref particles(bunch.get_local_particles());
    int kept = 0;
    int discarded = 0;
    int local_num = bunch.get_local_num();
    for (int part = 0; part < local_num; ++part) {
        bool try_discard = true;
        while (try_discard) {
            std::complex<double > u(particles[part][Bunch::x], particles[part][Bunch::y]);
            //std::cout << u.real() << "  " << u.imag() << std::endl;
            int index = 0;
            int size = vertices.size();
            double theta_sum = 0.0;
            while (index < size) {
                int index2 = index + 1;
                if (size == index2) index2 = 0;
                std::complex<double > v(vertices[index]);
                std::complex<double > w(vertices[index2]);
                double theta = arg((w - u) * conj(v - u));
                theta_sum += theta;
                ++index;
            }
            if (theta_sum / (2.0 * mconstants::pi) < 1.0e-12) {
                ++discarded;
                --local_num;
                std::cout << "polygon aperture: " << u.real() << "  " << u.imag() << std::endl;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else {
                ++kept;
                try_discard = false;
                //std::cout << u.real() << "  " << u.imag() << std::endl;
            }
        }
    }
    //std::cout << "kept = " << kept << ", discarded = " << discarded
    //        << std::endl;
    bunch.set_local_num(local_num);
}

Polygon_aperture_operation::~Polygon_aperture_operation()
{
}

const char Wire_elliptical_aperture_operation::type_name[] = "wire_elliptical_aperture";
const char Wire_elliptical_aperture_operation::attribute_name[] = "wire_elliptical";

Wire_elliptical_aperture_operation::Wire_elliptical_aperture_operation(
        Lattice_element const& element) :
    Aperture_operation(element)
{
    if (element.has_double_attribute("wire_elliptical_aperture_horizontal_radius")) {
        horizontal_radius = element.get_double_attribute(
                "wire_elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_horizontal_radius attribute");
    }
    if (element.has_double_attribute("wire_elliptical_aperture_vertical_radius")) {
        vertical_radius = element.get_double_attribute(
                "wire_elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_vertical_radius attribute");
    }
    if (element.has_double_attribute("wire_elliptical_aperture_wire_x")) {
         wire_x= element.get_double_attribute(
                "wire_elliptical_aperture_wire_x");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_x attribute");
    }
    if (element.has_double_attribute("wire_elliptical_aperture_wire_width")) {
        wire_width = element.get_double_attribute(
                "wire_elliptical_aperture_wire_width");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_width attribute");
    }
    if (element.has_double_attribute("wire_elliptical_aperture_gap")) {
        gap = element.get_double_attribute(
                "wire_elliptical_aperture_gap");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_gap attribute");
    }
}

const char *
Wire_elliptical_aperture_operation::get_type_name() const
{
    return type_name;
}

bool
Wire_elliptical_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (type_name == aperture_operation.get_type_name()) {
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
            && (wire_x
                    == wire_elliptical_aperture_operation.wire_x)
            && (wire_width
                    == wire_elliptical_aperture_operation.wire_width)
            && (gap
                    == wire_elliptical_aperture_operation.gap));
}

void
Wire_elliptical_aperture_operation::apply(Bunch & bunch)
{
    double h2 = horizontal_radius * horizontal_radius;
    double v2 = vertical_radius * vertical_radius;
    MArray2d_ref particles(bunch.get_local_particles());
    int kept = 0;
    int discarded = 0;
    int local_num = bunch.get_local_num();
    for (int part = 0; part < local_num; ++part) {
        bool try_discard = true;
        while (try_discard) {
            double x = particles[part][Bunch::x];
            double y = particles[part][Bunch::y];
            double scaled_r2 = x * x / h2 + y * y / v2;
            if (scaled_r2 > 1.0) {
                std::cout << "elliptical aperture: " << x << "  " << y << std::endl;
                ++discarded;
                --local_num;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else if ((x >= wire_x) && (x <= wire_x + wire_width)) {
                std::cout << "wire aperture: " << x << "  " << y << std::endl;
                ++discarded;
                --local_num;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else if (x >= wire_x + wire_width + gap) {
                std::cout << "gap aperture: " << x << "  " << y << std::endl;
                ++discarded;
                --local_num;
                if (part == local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = local_num;
                    particles[part][0] = particles[last][0];
                    particles[part][1] = particles[last][1];
                    particles[part][2] = particles[last][2];
                    particles[part][3] = particles[last][3];
                    particles[part][4] = particles[last][4];
                    particles[part][5] = particles[last][5];
                    particles[part][6] = particles[last][6];
                }
            } else {
                ++kept;
                try_discard = false;
            }
        }
    }
    //std::cout << "kept = " << kept << ", discarded = " << discarded
    //        << std::endl;
    bunch.set_local_num(local_num);
}

Wire_elliptical_aperture_operation::~Wire_elliptical_aperture_operation()
{
}
