#ifndef STEP_H_
#define STEP_H_

#include <list>
#include <boost/shared_ptr.hpp>
#include "synergia/utils/logger.h"
#include "synergia/simulation/operator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_train.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/utils/serialization.h"

struct Bunch_means
{
    double x_mean;
    double y_mean;
    double z_mean;
    double realnum;
    int bucket_index;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

class Step
{
private:
    Operators operators;
    std::list<double > time_fractions;
    double length;
    /// a list which contains informations about prevoius turns, the last element is the last (the earliest) turn stored
    /// it is updated at every step where impedance kick is applied
    /// the element is a vector of size num_bunches. The vector elements correspund to different bunches
    std::list< std::vector<Bunch_means> > stored_vbunches;

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
            Diagnosticss const& per_operation_diagnostics, Logger & logger);
    virtual void
    apply(Bunch_train & bunch_train, int verbosity,
            Train_diagnosticss const& per_operator_train_diagnosticss,
            Train_diagnosticss const& per_operation_train_diagnosticss,
            Logger & logger);
    Operators const&
    get_operators() const;
    Operators &
    get_operators();
    std::list<double> const&
    get_time_fractions() const;
    double
    get_length() const;
    std::list< std::vector<Bunch_means> >  const& get_stored_vbunches() const;
    virtual void
    print(int index) const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

typedef boost::shared_ptr<Step > Step_sptr;
typedef std::list<Step_sptr > Steps;

#endif /* STEP_H_ */
