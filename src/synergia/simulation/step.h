#ifndef STEP_H_
#define STEP_H_

#include "synergia/utils/logger.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/utils/cereal.h"


class Step
{

private:

    std::vector<std::pair<std::unique_ptr<Operator>, double>> operators;

    double length;
    std::vector<double> step_betas;

public:

    explicit Step(double length);

    template<typename OP>
    void append(OP const & opr, double time_fraction)
    { 
        operators.push_back(
                std::make_pair(std::make_unique<OP>(opr), time_fraction)); 
    }

    template<typename OP>
    void append(std::vector<OP> const & oprs, double time_fraction)
    { 
        for (auto const & opr : oprs) operators.push_back( 
            std::make_pair(std::make_unique<OP>(opr), time_fraction) ); 
    }

    void apply(Bunch_simulator & simulator, Logger & logger) const;
   
    double get_length() const;

#if 0
    void set_betas(double betax, double betay);
    std::vector<double > get_betas();

    void print(int index) const;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif
};

#endif /* STEP_H_ */
