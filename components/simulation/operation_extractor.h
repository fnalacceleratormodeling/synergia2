#ifndef OPERATION_EXTRACTOR_H_
#define OPERATION_EXTRACTOR_H_

#include "components/simulation/independent_operation.h"

class Operation_extractor
{
private:
    Chef_lattice_sptr chef_lattice_sptr;
    int map_order;
public:
    Operation_extractor(Chef_lattice_sptr const& chef_lattice_sptr, int map_order);
    Chef_lattice_sptr &
    get_chef_lattice_sptr();
    int
    get_map_order() const;
    virtual Independent_operations
    extract(Reference_particle const& reference_particle,
            Lattice_element_slices const& slices) = 0;
    virtual
    ~Operation_extractor();
};

typedef boost::shared_ptr<Operation_extractor > Operation_extractor_sptr;

class Chef_map_operation_extractor : public Operation_extractor
{
public:
    Chef_map_operation_extractor(Chef_lattice_sptr const& chef_lattice_sptr,
            int map_order);
    virtual Independent_operations
    extract(Reference_particle const& reference_particle,
            Lattice_element_slices const& slices);
};

class Chef_propagate_operation_extractor : public Operation_extractor
{
public:
    Chef_propagate_operation_extractor(Chef_lattice_sptr const& chef_lattice_sptr,
            int map_order);
    virtual Independent_operations
    extract(Reference_particle const& reference_particle,
            Lattice_element_slices const& slices);
};

class Mixed_chef_operation_extractor : public Operation_extractor
{
public:
    Mixed_chef_operation_extractor(Chef_lattice_sptr const& chef_lattice_sptr,
            int map_order);
    virtual Independent_operations
    extract(Reference_particle const& reference_particle,
            Lattice_element_slices const& slices);
};

class Operation_extractor_map
{
private:
    std::map<std::string, Operation_extractor_sptr > extractor_map;
public:
    Operation_extractor_map();
    void
    set_extractor(std::string const& name,
            Operation_extractor_sptr const& operation_extractor);
    Operation_extractor_sptr
    get_extractor(std::string const& name);
    ~Operation_extractor_map();
};

#endif /* OPERATION_EXTRACTOR_H_ */
