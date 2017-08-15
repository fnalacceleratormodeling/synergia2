#ifndef STEP_H_
#define STEP_H_

#include <list>
#include <boost/shared_ptr.hpp>
#include "synergia/utils/logger.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/propagate_actions.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_train.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/utils/serialization.h"



class Step
{
private:
    Operators operators;
    std::list<double > time_fractions;
    double length;
    std::vector<double > step_betas;
   

public:
    Step(double length);
    // Default constructor for serialization use only
    Step();
    void
    append(Operator_sptr operator_sptr, double time_fraction);
    void
    append(Operators const& operators, double time_fraction);
    virtual void
    apply(Bunch & bunch, int verbosity,
            Diagnosticss const& per_operator_diagnostics,
            Diagnosticss const& per_operation_diagnostics, Stepper & stepper, Logger & logger);

    virtual void
    apply(Bunch_train & bunch_train, int verbosity,
        Train_diagnosticss const& per_operator_train_diagnosticss,
        Train_diagnosticss const& per_operation_train_diagnosticss, 
         Propagate_actions * propagate_actions_ptr, Stepper & stepper, int step_count,  int turn,  
         Logger & logger);
    
    
    Operators const&
    get_operators() const;
    Operators &
    get_operators();
    std::list<double> const&
    get_time_fractions() const;
    double
    get_length() const;
    void
    set_betas(double betax, double betay);
    std::vector<double >
    get_betas();
    virtual void
    print(int index) const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

typedef boost::shared_ptr<Step > Step_sptr; // syndoc:include
typedef std::list<Step_sptr > Steps; // syndoc:include

#endif /* STEP_H_ */
