#ifndef INDEPENDENT_OPERATOR_H_
#define INDEPENDENT_OPERATOR_H_

#include "synergia/simulation/operator.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/simulation/independent_operation.h"
#include "synergia/simulation/operation_extractor.h"
#include "synergia/simulation/aperture_operation_extractor.h"
#include "synergia/foundation/multi_diagnostics.h"
#
class Independent_operator : public Operator
{
private:
    Lattice_element_slices slices;
    Independent_operations operations;
    std::list<long int > operations_revisions;
    Reference_particle operations_reference_particle;
    Operation_extractor_map_sptr operation_extractor_map_sptr;
    Aperture_operation_extractor_map_sptr aperture_operation_extractor_map_sptr;
    bool have_operations;
public:
    Independent_operator(
            std::string const& name,
            Operation_extractor_map_sptr operation_extractor_map_sptr,
            Aperture_operation_extractor_map_sptr aperture_operation_extractor_map_sptr);
    /// Default constructor for serialization use only
    Independent_operator();
    void
    append_slice(Lattice_element_slice_sptr slice_sptr);
    Lattice_element_slices const&
    get_slices() const;
    void
    update_operations(Reference_particle const& reference_particle);
    bool
    need_update(Reference_particle const& reference_particle, int verbosity,
            Logger & logger);
    Independent_operations const&
    get_operations() const;
    Independent_operations &
    get_operations();
    using Operator::apply;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Diagnosticss const& per_operation_diagnosticss, Logger & logger);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger, Multi_diagnostics & diagnostics);
    virtual void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    static const char type_name[];
    virtual
    ~Independent_operator();
};
BOOST_CLASS_EXPORT_KEY(Independent_operator);
typedef boost::shared_ptr<Independent_operator > Independent_operator_sptr; // syndoc:include
typedef std::list<Independent_operator_sptr > Independent_operators; // syndoc:include



#endif
