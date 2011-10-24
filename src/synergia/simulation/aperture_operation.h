#ifndef APERTURE_OPERATION_H_
#define APERTURE_OPERATION_H_

#include "synergia/simulation/independent_operation.h"

class Aperture_operation : public Independent_operation
{
public:
    Aperture_operation(Lattice_element const& element);
    virtual const char *
    get_type_name() const = 0;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const = 0;
    virtual void
    apply(Bunch & bunch)=0;
    virtual
    ~Aperture_operation();
};

typedef boost::shared_ptr<Aperture_operation > Aperture_operation_sptr;

/// A circular aperture with radius in meters determined by the
/// Lattice_element attribute "circular_aperture_radius".
/// If the radius is not defined, the default value of 1000.0 m will
/// be used.
class Circular_aperture_operation : public Aperture_operation
{
private:
    double radius;
public:
    static const double default_radius;
    static const char type_name[];
    static const char attribute_name[];
    Circular_aperture_operation(Lattice_element const& element);
    virtual const char *
    get_type_name() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    virtual bool
    operator==(Circular_aperture_operation const& circular_aperture_operation) const;
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Circular_aperture_operation();
};

/// An elliptical aperture with horizontal and vertical radii in meters
/// determined by the Lattice_element_attributes
/// "elliptical_aperture_horizontal_radius" and
/// "elliptical_aperture_vertical_radius", respectively.
/// Both radii must be specified. Failing to do so will cause an
/// exception.
class Elliptical_aperture_operation : public Aperture_operation
{
private:
    double horizontal_radius, vertical_radius;
public:
    static const char type_name[];
    static const char attribute_name[];
    Elliptical_aperture_operation(Lattice_element const& element);
    virtual const char *
    get_type_name() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Elliptical_aperture_operation();
};

/// A rectangular aperture with horizontal and vertical dimensions in meters
/// determined by the Lattice_element_attributes
/// "rectangular_aperture_width" and
/// "rectangular_aperture_height", respectively.
/// Both dimensions must be specified. Failing to do so will cause an
/// exception.
class Rectangular_aperture_operation : public Aperture_operation
{
private:
    double width, height;
public:
    static const char type_name[];
    static const char attribute_name[];
    Rectangular_aperture_operation(Lattice_element const& element);
    virtual const char *
    get_type_name() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Rectangular_aperture_operation();
};

#endif /* APERTURE_OPERATION_H_ */
