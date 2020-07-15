#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>
#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_train.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/logger.h"
#include "synergia/simulation/propagate_actions.h"

class Step;

class Operator
{
private:
    std::string name, type;
public:
    Operator(std::string const& name, std::string const& type);
    /// Default constructor for serialization use only
    Operator();
    std::string const&
    get_name() const;
    std::string const&
    get_type() const;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Diagnosticss const& per_operation_diagnosticss, Logger & logger) = 0;
    virtual void
    apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
            Train_diagnosticss const& per_operation_train_diagnosticss, Logger & logger);
            
    virtual void
    apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
            Train_diagnosticss const& per_operation_train_diagnosticss, 
            Propagate_actions * propagate_actions_ptr, Stepper & stepper, int step_count,  int turn, 
            Logger & logger);
                    
    virtual void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Operator();
};
typedef boost::shared_ptr<Operator > Operator_sptr; // syndoc:include
typedef std::list<Operator_sptr > Operators; // syndoc:include

#endif /* OPERATOR_H_ */
