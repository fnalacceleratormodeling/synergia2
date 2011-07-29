#ifndef STEP_H_
#define STEP_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/simulation/operator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/multi_diagnostics.h"

struct Bunch_means
{
double x_mean;
double y_mean;
double z_mean;
double n_part;
};


class Step
{
private:
    Operators operators;
    std::list<double > time_fractions;
    double length;
    std::list<Bunch_means> stored_bunches;
public:
    Step(double length);
    void
    append(Operator_sptr operator_sptr, double time_fraction);
    void
    append(Operators const& operators, double time_fraction);
    virtual void
    apply(Bunch & bunch);
    virtual void
    apply(Bunch & bunch, Multi_diagnostics & diagnostics);
    Operators const&
    get_operators() const;
    std::list<double> const&
    get_time_fractions() const;
    double
    get_length() const;
    std::list<Bunch_means>  get_stored_bunches() const;
    virtual void
    print(int index) const;
};

typedef boost::shared_ptr<Step > Step_sptr;
typedef std::list<Step_sptr > Steps;

#endif /* STEP_H_ */
