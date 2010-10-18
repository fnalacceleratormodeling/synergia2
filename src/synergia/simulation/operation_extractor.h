#ifndef OPERATION_EXTRACTOR_H_
#define OPERATION_EXTRACTOR_H_

#include "synergia/simulation/independent_operation.h"

class Operation_extractor
{
private:
    Chef_lattice_sptr chef_lattice_sptr;
    int map_order;
public:
    Operation_extractor(Chef_lattice_sptr chef_lattice_sptr, int map_order);
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
    Chef_map_operation_extractor(Chef_lattice_sptr chef_lattice_sptr,
            int map_order);
    virtual Independent_operations
    extract(Reference_particle const& reference_particle,
            Lattice_element_slices const& slices);
};

class Chef_propagate_operation_extractor : public Operation_extractor
{
public:
    Chef_propagate_operation_extractor(Chef_lattice_sptr chef_lattice_sptr,
            int map_order);
    virtual Independent_operations
    extract(Reference_particle const& reference_particle,
            Lattice_element_slices const& slices);
};

class Chef_mixed_operation_extractor : public Operation_extractor
{
public:
    Chef_mixed_operation_extractor(Chef_lattice_sptr chef_lattice_sptr,
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
            Operation_extractor_sptr operation_extractor);
    Operation_extractor_sptr &
    get_extractor(std::string const& name);
    std::list<std::string >
    get_extractor_names() const;
    ~Operation_extractor_map();
};

typedef boost::shared_ptr<Operation_extractor_map > Operation_extractor_map_sptr;

const char default_operation_extractor_name[] = "default";
const char chef_map_operation_extractor_name[] = "chef_map";
const char chef_propagate_operation_extractor_name[] = "chef_propagate";
const char chef_mixed_operation_extractor_name[] = "chef_mixed";

#endif /* OPERATION_EXTRACTOR_H_ */
