#ifndef APERTURE_OPERATION_EXTRACTOR_H_
#define APERTURE_OPERATION_EXTRACTOR_H_

#include "synergia/simulation/aperture_operation.h"

class Aperture_operation_extractor
{
public:
    Aperture_operation_extractor();
    virtual Aperture_operation_sptr
    extract(Lattice_element const& element) = 0;
    virtual
    ~Aperture_operation_extractor();
};

typedef boost::shared_ptr<Aperture_operation_extractor >
        Aperture_operation_extractor_sptr;

template<typename T>
    class Generic_aperture_extractor : public Aperture_operation_extractor
    {
    public:
        Generic_aperture_extractor()
        {
        }
        virtual Aperture_operation_sptr
        extract(Lattice_element const& element)
        {
            return boost::shared_ptr<T >(new T(element));
        }
        virtual
        ~Generic_aperture_extractor()
        {
        }
    };

class Aperture_operation_extractor_map
{
private:
    std::map<std::string, Aperture_operation_extractor_sptr > extractor_map;
public:
    Aperture_operation_extractor_map();
    void
    set_extractor(std::string const& name,
            Aperture_operation_extractor_sptr operation_extractor);
    Aperture_operation_extractor_sptr
    get_extractor(std::string const& name);
    std::list<std::string >
    get_extractor_names() const;
    ~Aperture_operation_extractor_map();
};

typedef boost::shared_ptr<Aperture_operation_extractor_map >
        Aperture_operation_extractor_map_sptr;

#endif /* APERTURE_OPERATION_EXTRACTOR_H_ */
