#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/lattice/lattice_element_slice.h"

#include "synergia/utils/cereal.h"
#include "synergia/utils/logger.h"

class Step;

class Operator
{

private:

    std::string name;
    std::string type;

private:

    virtual void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) = 0;

public:

    Operator() : name(), type() { }

    Operator(std::string const & name, std::string const & type)
        : name(name), type(type)
    { }

    virtual ~Operator() = default;

    std::string const & get_name() const { return name; }
    std::string const & get_type() const { return type; }

    void apply(Bunch_simulator & simulator, double time_step, Logger & logger)
    { apply_impl(simulator, time_step, logger); }
                    
    //virtual void print() const;

#if 0
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    { ar & CEREAL_NVP(name) & CEREAL_NVP(type); }
#endif
};


class Collective_operator : public Operator
{
private:

public:

    Collective_operator() : Operator() { }

    Collective_operator(std::string const & name)
        : Operator(name, "collective")
    { }

#if 0
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    { }
#endif
};

class Dummy_collective_operator : public Collective_operator
{
private:

    virtual void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) override
    { }

public:

    Dummy_collective_operator() : Collective_operator() { }

    Dummy_collective_operator(std::string const& name)
        : Collective_operator("dummy_collective")
    { }

#if 0
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    { }
#endif
};


class Independent_operator : public Operator
{

private:

    Independent_operator();

private:

    std::vector<Lattice_element_slice> slices;
    //std::vector<Independent_operations> operations;

private:

#if 0
    Lattice_element_slices slices;
    Independent_operations operations;
    std::list<long int > operations_revisions;
    Reference_particle operations_reference_particle;
    Operation_extractor_map_sptr operation_extractor_map_sptr;
    Aperture_operation_extractor_map_sptr aperture_operation_extractor_map_sptr;
    bool have_operations;
#endif

public:

    Independent_operator(std::string const & name);

    void append_slice(Lattice_element_slice const & slice);
    std::vector<Lattice_element_slice> const & get_slices() const;

#if 0
    void update_operations(
            Reference_particle const& reference_particle);

    bool need_update(
            Reference_particle const & reference_particle, 
            int verbosity, Logger & logger);

    Independent_operations const & get_operations() const;
    Independent_operations       & get_operations();
#endif

    void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) override;

#if 0
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif
};

#endif /* OPERATOR_H_ */
