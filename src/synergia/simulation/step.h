#ifndef STEP_H_
#define STEP_H_

#include "synergia/utils/logger.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/bunch_simulator.h"


class Step
{

private:

    std::vector<std::unique_ptr<Operator>> operators;

    double length;
    std::vector<double> step_betas;

public:

    explicit Step(double length);

    template<class... Args>
    Independent_operator & append_independent(Args &&... args)
    { 
        operators.emplace_back(
                std::make_unique<Independent_operator>(std::forward<Args>(args)...) );
        return *(dynamic_cast<Independent_operator*>(operators.back().get()));
    }

    void append_collective(std::unique_ptr<CO_options> const & co_ops)
    {
        operators.emplace_back(co_ops->create_operator());
    }

    void apply(Bunch_simulator & simulator, Logger & logger) const;
    void create_operations(Lattice const & lattice);
   
    double get_length() const
    { return length; }

    void print(Logger & logger) const
    { 
        logger(LoggerV::DEBUG) << "step length = " << length << "\n";
        for(auto const & opr : operators) opr->print(logger); 
        logger(LoggerV::DEBUG) << "\n";
    }

#if 0
    void set_betas(double betax, double betay);
    std::vector<double > get_betas();

    void print(int index) const;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif
};

#endif /* STEP_H_ */
