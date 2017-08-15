#ifndef APERTURE_OPERATION_H_
#define APERTURE_OPERATION_H_

#include "synergia/simulation/independent_operation.h"

class Aperture_operation : public Independent_operation
{
private:
    Lattice_element_slice_sptr slice_sptr;
    double x_offset, y_offset;
public:
    static const char charge_attribute[];
    Aperture_operation(Lattice_element_slice_sptr slice_sptr);
    // Default constructor for serialization use only
    Aperture_operation();
    virtual const char *
    get_aperture_type() const = 0;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const = 0;
    template<typename T>
        void
        apply_impl(T & t, Bunch & bunch, int verbosity, Logger & logger);
    template<typename T>
        void
        dump_particles(T & t, Bunch & bunch, int verbosity, Logger & logger);
    inline double
    get_x_offset() const
    {
        return x_offset;
    }
    inline double
    get_y_offset() const
    {
        return y_offset;
    }
    Lattice_element_slice_sptr 
    get_slice_sptr() const;
    virtual void
    apply(Bunch & bunch, int verbosity, Logger & logger)=0;
    void
    deposit_charge(double charge);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Aperture_operation)

typedef boost::shared_ptr<Aperture_operation > Aperture_operation_sptr; // syndoc:include
typedef std::list<Aperture_operation_sptr > Aperture_operation_sptrs; // syndoc:include

/// An aperture to remove all particles with infinite and/or NaN coordinates.
class Finite_aperture_operation : public Aperture_operation
{
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Finite_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    // Default constructor for serialization use only
    Finite_aperture_operation();
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
    virtual void
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Finite_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Finite_aperture_operation)

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
    // Default constructor for serialization use only
    Circular_aperture_operation();
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
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Circular_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Circular_aperture_operation)

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
    // Default constructor for serialization use only
    Elliptical_aperture_operation();
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
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Elliptical_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Elliptical_aperture_operation)
;

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
    // Default constructor for serialization use only
    Rectangular_aperture_operation();
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
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Rectangular_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Rectangular_aperture_operation)
;

/// A rectangular aperture with quarter-circular ears on top of a
/// rectangular strip. The rectangular
/// horizontal and vertical dimensions in meters
/// determined by the Lattice_element_attributes
/// "rectangular_aperture_width" and
/// "rectangular_aperture_height", respectively.  The quarter circle radius
/// is given by "ear_radius" and the vertical location of the
/// the center of the ear is given by option "rectangular_aperture_ear_offset".  The aperture
/// has four-fold symmetry.
/// All dimensions must be specified. Failing to do so will cause an
/// exception.
class Rectangular_with_ears_aperture_operation : public Aperture_operation
{
private:
    double width, height, ear_offset, radius;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Rectangular_with_ears_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    // Default constructor for serialization use only
    Rectangular_with_ears_aperture_operation();
    virtual const char *
            get_aperture_type() const;
    virtual bool
            operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Rectangular_with_ears_aperture_operation const& rectangular_with_ears_aperture_operation) const;
    bool
            operator()(MArray2d_ref & particles, int part);
    virtual void
            apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
    void
            serialize(Archive & ar, const unsigned int version);
    virtual
            ~Rectangular_with_ears_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Rectangular_with_ears_aperture_operation)
;

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
    int num_vertices;
    std::vector<std::complex<double > > vertices;
    double min_radius2;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Polygon_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    // Default constructor for serialization use only
    Polygon_aperture_operation();
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
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Polygon_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Polygon_aperture_operation)

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
    // Default constructor for serialization use only
    Wire_elliptical_aperture_operation();
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
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Wire_elliptical_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Wire_elliptical_aperture_operation)

/// A Lambertson aperture with radius in meters determined by the
/// Lattice_element attribute "lambertson_aperture_radius".
/// If the radius is not defined, the default value of 1000.0 m will
/// be used.
class Lambertson_aperture_operation : public Aperture_operation
{
private:
    double radius;
public:
    static const char aperture_type[];
    static const char attribute_name[];
    Lambertson_aperture_operation(Lattice_element_slice_sptr slice_sptr);
    // Default constructor for serialization use only
    Lambertson_aperture_operation();
    virtual const char *
    get_aperture_type() const;
    virtual bool
    operator==(Aperture_operation const& aperture_operation) const;
    bool
            operator==(
                    Lambertson_aperture_operation const& lambertson_aperture_operation) const;
    bool
    operator()(MArray2d_ref & particles, int part);
    virtual void
    apply(Bunch & bunch, int verbosit, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Lambertson_aperture_operation();
};
BOOST_CLASS_EXPORT_KEY(Lambertson_aperture_operation)


#include "synergia/simulation/aperture_operation.tcc"

#endif /* APERTURE_OPERATION_H_ */
