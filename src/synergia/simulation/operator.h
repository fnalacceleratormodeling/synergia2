#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>

#include "synergia/lattice/lattice_element_slice.h"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_operation.h"
#include "synergia/simulation/operation_extractor.h"

#include "synergia/utils/logger.h"

class Step;

class Operator
{

private:

    std::string name;
    std::string type;

    double time_fraction;

private:

    virtual void create_operations_impl(
            Lattice const & lattice) = 0;

    virtual void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) = 0;

public:

    Operator() : name(), type(), time_fraction(1.0) { }

    Operator(std::string const & name, std::string const & type, double time)
        : name(name), type(type), time_fraction(time)
    { }

    virtual ~Operator() = default;
    virtual Operator * clone() const = 0;

    std::string const & get_name() const { return name; }
    std::string const & get_type() const { return type; }
    double get_time_fraction() const { return time_fraction; }

    void create_operations(Lattice const & lattice)
    { create_operations_impl(lattice); }

    void apply(Bunch_simulator & simulator, double time_step, Logger & logger)
    { apply_impl(simulator, time_step * time_fraction, logger); }
                    
    //virtual void print() const;
};

class Collective_operator : public Operator
{
private:

    void create_operations_impl(
            Lattice const & lattice) final
    { }

public:

    Collective_operator() : Operator() { }

    Collective_operator(std::string const & name, double time)
        : Operator(name, "collective", time)
    { }
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

    Dummy_collective_operator(std::string const& name, double time)
        : Collective_operator("dummy_collective", time)
    { }

    Dummy_collective_operator * clone() const override
    { return new Dummy_collective_operator(*this); }
};


class Independent_operator : public Operator
{
private:

    void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) override;

    void create_operations_impl(
            Lattice const & lattice) override;

    bool need_update(
            Reference_particle const & ref, 
            Logger & logger);

    void update_operations(
            Reference_particle const & ref);

private:

    std::vector<Lattice_element_slice> slices;
    std::vector<std::unique_ptr<Independent_operation>> operations;

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

    Independent_operator(std::string const & name, double time);

    Independent_operator * clone() const override
    { return nullptr; }

    template<class... Args>
    Independent_operator & append_slice(Args && ... args)
    { slices.emplace_back(std::forward<Args>(args)...); return *this; }

    std::vector<Lattice_element_slice> const & get_slices() const
    { return slices; }
};

#endif /* OPERATOR_H_ */


