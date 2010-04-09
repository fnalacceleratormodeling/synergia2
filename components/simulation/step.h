#ifndef STEP_H_
#define STEP_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "components/simulation/operator.h"
#include "components/bunch/bunch.h"

class Step
{
private:
    Operators operators;
public:
    Step();
    void
    append(Operator_sptr operator_sptr);
    void
    append(Operators const& operators);
    void
    apply(Bunch & bunch);
    Operators const&
    get_operators() const;
    void
    print(int index) const;
};

typedef boost::shared_ptr<Step > Step_sptr;
typedef std::list<Step_sptr > Steps;

#endif /* STEP_H_ */
