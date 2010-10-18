#ifndef STEP_H_
#define STEP_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/simulation/operator.h"
#include "synergia/bunch/bunch.h"

class Step
{
private:
    Operators operators;
    std::list<double > time_fractions;
public:
    Step();
    void
    append(Operator_sptr operator_sptr, double time_fraction);
    void
    append(Operators const& operators, double time_fraction);
    virtual void
    apply(Bunch & bunch);
    Operators const&
    get_operators() const;
    std::list<double> const&
    get_time_fractions() const;
    virtual void
    print(int index) const;
};

typedef boost::shared_ptr<Step > Step_sptr;
typedef std::list<Step_sptr > Steps;

#endif /* STEP_H_ */
