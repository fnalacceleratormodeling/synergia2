#ifndef APERTURE_OPERATION_H_
#define APERTURE_OPERATION_H_

#include "synergia/simulation/independent_operation.h"

const double default_circular_aperture_radius = 1000.0;
const char circular_aperture_type_name[] = "circular_aperture";
const char circular_aperture_attribute_name[] = "circular";
/// A circular aperture with radius in meters determined by the
/// Lattice_element attribute "circular_aperture_radius".
/// If the radius is not defined, the default value of 1000.0 m will
/// be used.
class Circular_aperture_operation : public Independent_operation
{
private:
    double radius;
public:
    Circular_aperture_operation(Lattice_element const& element);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Circular_aperture_operation();
};

const char elliptical_aperture_type_name[] = "elliptical_aperture";
const char elliptical_aperture_attribute_name[] = "elliptical";
/// An elliptical aperture with horizontal and vertical radii in meters
/// determined by the Lattice_element_attributes
/// "elliptical_aperture_horizontal_radius" and
/// "elliptical_aperture_vertical_radius", respectively.
/// Both radii must be specified. Failing to do so will cause an
/// exception.
class Elliptical_aperture_operation : public Independent_operation
{
private:
    double horizontal_radius, vertical_radius;
public:
    Elliptical_aperture_operation(Lattice_element const& element);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Elliptical_aperture_operation();
};

const char rectangular_aperture_type_name[] = "rectangular_aperture";
const char rectangular_aperture_attribute_name[] = "rectangular";
/// A rectangular aperture with horizontal and vertical dimensions in meters
/// determined by the Lattice_element_attributes
/// "rectangular_aperture_width" and
/// "rectangular_aperture_height", respectively.
/// Both dimensions must be specified. Failing to do so will cause an
/// exception.
class Rectangular_aperture_operation : public Independent_operation
{
private:
    double width, height;
public:
    Rectangular_aperture_operation(Lattice_element const& element);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Rectangular_aperture_operation();
};

#endif /* APERTURE_OPERATION_H_ */
