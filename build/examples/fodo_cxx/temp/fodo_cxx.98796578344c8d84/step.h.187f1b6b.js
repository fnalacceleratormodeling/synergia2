var data = {lines:[
{"lineNum":"    1","line":"#ifndef STEP_H_"},
{"lineNum":"    2","line":"#define STEP_H_"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":"#include \"synergia/utils/logger.h\""},
{"lineNum":"    5","line":"#include \"synergia/simulation/operator.h\""},
{"lineNum":"    6","line":"#include \"synergia/simulation/collective_operator_options.h\""},
{"lineNum":"    7","line":"#include \"synergia/simulation/bunch_simulator.h\""},
{"lineNum":"    8","line":""},
{"lineNum":"    9","line":"class Propagator;"},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"class Step","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   12","line":"{"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"private:"},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"    std::vector<std::shared_ptr<Operator>> operators;"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"    double length;"},
{"lineNum":"   19","line":"    std::vector<double> step_betas;"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"public:"},
{"lineNum":"   22","line":""},
{"lineNum":"   23","line":"    explicit Step(double length);"},
{"lineNum":"   24","line":""},
{"lineNum":"   25","line":"    template<class... Args>"},
{"lineNum":"   26","line":"    Independent_operator & append_independent(Args &&... args)"},
{"lineNum":"   27","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   28","line":"        operators.emplace_back(","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   29","line":"                std::make_shared<Independent_operator>(std::forward<Args>(args)...) );"},
{"lineNum":"   30","line":"        return *(dynamic_cast<Independent_operator*>(operators.back().get()));","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   31","line":"    }"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"    void append_collective(std::unique_ptr<CO_options> const & co_ops)"},
{"lineNum":"   34","line":"    { operators.emplace_back(co_ops->create_operator()); }"},
{"lineNum":"   35","line":""},
{"lineNum":"   36","line":"    void append(std::shared_ptr<Operator> col_op)"},
{"lineNum":"   37","line":"    { operators.push_back(col_op); }","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"    void apply(Bunch_simulator & simulator, Logger & logger) const;"},
{"lineNum":"   40","line":"    void create_operations(Lattice const & lattice);"},
{"lineNum":"   41","line":""},
{"lineNum":"   42","line":"    double get_length() const"},
{"lineNum":"   43","line":"    { return length; }"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"    void print(Logger & logger) const"},
{"lineNum":"   46","line":"    {"},
{"lineNum":"   47","line":"        logger(LoggerV::DEBUG) << \"step length = \" << length << \"\\n\";"},
{"lineNum":"   48","line":"        for(auto const & opr : operators) opr->print(logger);"},
{"lineNum":"   49","line":"        logger(LoggerV::DEBUG) << \"\\n\";"},
{"lineNum":"   50","line":"    }"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"#if 0"},
{"lineNum":"   53","line":"    void set_betas(double betax, double betay);"},
{"lineNum":"   54","line":"    std::vector<double > get_betas();"},
{"lineNum":"   55","line":""},
{"lineNum":"   56","line":"    void print(int index) const;"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"    template<class Archive>"},
{"lineNum":"   59","line":"    void serialize(Archive & ar, const unsigned int version);"},
{"lineNum":"   60","line":"#endif"},
{"lineNum":"   61","line":""},
{"lineNum":"   62","line":"    friend class Propagator;"},
{"lineNum":"   63","line":"};"},
{"lineNum":"   64","line":""},
{"lineNum":"   65","line":"#endif /* STEP_H_ */"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 5, "covered" : 0,};
var merged_data = [];
