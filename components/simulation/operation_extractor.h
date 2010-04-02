#ifndef OPERATION_EXTRACTOR_H_
#define OPERATION_EXTRACTOR_H_

#include "components/simulation/independent_operation.h"

class Operation_extractor
{
public:
    virtual Independent_operations
            extract(Lattice_element_slices const& slices,
                    Chef_lattice & chef_lattice) = 0;
    virtual
    ~Operation_extractor()
    {
    }
    ;
};

typedef boost::shared_ptr<Operation_extractor > Operation_extractor_sptr;
typedef std::map<std::string, Operation_extractor_sptr > Operation_extractor_map;

class Mixed_chef_operation_extractor : public Operation_extractor
{
public:
    virtual Independent_operations
    extract(Lattice_element_slices const& slices, Chef_lattice & chef_lattice);
};

class Chef_propagate_operation_extractor : public Operation_extractor
{
public:
    virtual Independent_operations
    extract(Lattice_element_slices const& slices, Chef_lattice & chef_lattice);
};

class Chef_map_operation_extractor : public Operation_extractor
{
public:
    virtual Independent_operations
    extract(Lattice_element_slices const& slices, Chef_lattice & chef_lattice);
};

#endif /* OPERATION_EXTRACTOR_H_ */
