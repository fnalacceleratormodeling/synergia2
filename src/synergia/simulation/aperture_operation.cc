#include "aperture_operation.h"
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
    //        std::cout << "kept = " << kept << ", discarded = " << discarded
    //                << std::endl;
    bunch.set_local_num(local_num);
    bunch.update_total_num();
   // std::cout<<"circular aperture applied"<<std::endl;
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
    //        std::cout << "kept = " << kept << ", discarded = " << discarded
    //                << std::endl;
    bunch.set_local_num(local_num);
    bunch.update_total_num();
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
    //        std::cout << "kept = " << kept << ", discarded = " << discarded
    //                << std::endl;
    bunch.set_local_num(local_num);
    bunch.update_total_num();
  //  std::cout<<"rectangular aperture applied"<<std::endl;
}

Rectangular_aperture_operation::~Rectangular_aperture_operation()
{
}
