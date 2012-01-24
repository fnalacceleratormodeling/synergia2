#ifndef APERTURE_OPERATION_H_
#define APERTURE_OPERATION_H_

#include "synergia/simulation/independent_operation.h"

class Aperture_operation : public Independent_operation
{
private:
    Lattice_element_slice_sptr slice_sptr;
public:
    static const char charge_attribute[];
    Aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const = 0;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const = 0;
    template<typename T>
        void
        apply_impl(T& t, Bunch & bunch);
    virtual void
    apply(Bunch & bunch)=0;
    void
    deposit_charge(double charge);
    virtual
    ~Aperture_operation();
};

typedef boost::shared_ptr<Aperture_operation > Aperture_operation_sptr;
typedef std::list<Aperture_operation_sptr > Aperture_operation_sptrs;

/// An aperture to remove all particles with infinite and/or NaN coordinates.
class Finite_aperture_operation : public Aperture_operation
{
private:
    double radius;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Finite_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Finite_aperture_operation();
};

/// A circular aperture with radius in meters determined by the
/// Lattice_element attribute "circular_aperture_radius".
/// If the radius is not defined, the default value of 1000.0 m will
/// be used.
class Circular_aperture_operation : public Aperture_operation
{
private:
    double radius, radius2;
public:
    static const double default_radius;
    static const char aperture_type[];
    static const char attribute_name[];
    Circular_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Circular_aperture_operation const& circular_aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
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
    double h2, v2;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Elliptical_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Elliptical_aperture_operation const& Elliptical_aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
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
    static const char aperture_type[];
    static const char attribute_name[];
    Rectangular_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Rectangular_aperture_operation const& rectangular_aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Rectangular_aperture_operation();
};

/// A polygon aperture with vertices
/// determined by the Lattice_element_attributes
/// "pax1", "pay1", "pax2", "pay2", and so on.
/// And it also requires "the_number_of_vertices", which determines the number
/// of vertices and must be greter than and equal to 3.
/// Must have at least 3 vertcies. Failing to do so will cause an
/// exception.
class Polygon_aperture_operation : public Aperture_operation
{
private:
    int vertices_num;
    std::vector<std::complex<double > > vertices;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Polygon_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Polygon_aperture_operation const& polygon_aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Polygon_aperture_operation();
};

/// An wire_elliptical aperture with horizontal and vertical radii in meters
/// determined by the Lattice_element_attributes
/// "wire_elliptical_aperture_horizontal_radius" and
/// "wire_elliptical_aperture_vertical_radius", respectively.
/// Both radii must be specified.
/// Also needs Lattice_element_attributes
/// "wire_elliptical_aperture_wire_x",
/// "wire_elliptical_aperture_wire_width", and
/// "wire_elliptical_aperture_gap." Failing to do so will cause an
/// exception.
class Wire_elliptical_aperture_operation : public Aperture_operation
{
private:
    double horizontal_radius, vertical_radius;
    double h2, v2;
    double wire_x, wire_width, gap;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Wire_elliptical_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Wire_elliptical_aperture_operation const& Wire_elliptical_aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Wire_elliptical_aperture_operation();
};

#include "synergia/simulation/aperture_operation.tcc"

#endif /* APERTURE_OPERATION_H_ */
