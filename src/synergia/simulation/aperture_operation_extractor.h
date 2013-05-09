#ifndef APERTURE_OPERATION_EXTRACTOR_H_
#define APERTURE_OPERATION_EXTRACTOR_H_

#include "synergia/simulation/aperture_operation.h"

class Aperture_operation_extractor {
public:
    Aperture_operation_extractor();
    virtual Aperture_operation_sptr
    extract(Lattice_element_slice_sptr slice) = 0;
	virtual
	~Aperture_operation_extractor();
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
};

typedef boost::shared_ptr<Aperture_operation_extractor> Aperture_operation_extractor_sptr; // syndoc:include

template<typename T>
class Generic_aperture_extractor: public Aperture_operation_extractor {
public:
	Generic_aperture_extractor() {
	}
        virtual Aperture_operation_sptr
        extract(Lattice_element_slice_sptr slice)
        {
            return boost::shared_ptr<T >(new T(slice));
	}
	virtual ~Generic_aperture_extractor() {
	}
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation_extractor);
	}
};

typedef Generic_aperture_extractor<Circular_aperture_operation> Circular_extractor; // syndoc:include
BOOST_CLASS_EXPORT_KEY(Circular_extractor);

typedef Generic_aperture_extractor<Elliptical_aperture_operation> Elliptical_extractor; // syndoc:include
BOOST_CLASS_EXPORT_KEY(Elliptical_extractor);

typedef Generic_aperture_extractor<Rectangular_aperture_operation> Rectangular_extractor; // syndoc:include
BOOST_CLASS_EXPORT_KEY(Rectangular_extractor);

typedef Generic_aperture_extractor<Polygon_aperture_operation >
        Polygon_extractor; // syndoc:include
BOOST_CLASS_EXPORT_KEY(Polygon_extractor);

typedef Generic_aperture_extractor<Wire_elliptical_aperture_operation >
        Wire_elliptical_extractor; // syndoc:include
BOOST_CLASS_EXPORT_KEY(Wire_elliptical_extractor);

typedef Generic_aperture_extractor<Lambertson_aperture_operation >
        Lambertson_extractor; // syndoc:include
BOOST_CLASS_EXPORT_KEY(Lambertson_extractor);

class Aperture_operation_extractor_map {
private:
	std::map<std::string, Aperture_operation_extractor_sptr> extractor_map;
public:
	Aperture_operation_extractor_map();
	void
	set_extractor(std::string const& name,
			Aperture_operation_extractor_sptr operation_extractor);
	Aperture_operation_extractor_sptr
	get_extractor(std::string const& name);
	std::list<std::string>
	get_extractor_names() const;
	~Aperture_operation_extractor_map();
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
};

typedef boost::shared_ptr<Aperture_operation_extractor_map> Aperture_operation_extractor_map_sptr; // syndoc:include

#endif /* APERTURE_OPERATION_EXTRACTOR_H_ */

