var data = {lines:[
{"lineNum":"    1","line":"#ifndef SPLIT_OPERATOR_STEPPER_ELEMENTS_H_"},
{"lineNum":"    2","line":"#define SPLIT_OPERATOR_STEPPER_ELEMENTS_H_"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":"#include \"synergia/simulation/stepper.h\""},
{"lineNum":"    5","line":"#include \"synergia/collective/dummy_collective_operator.h\""},
{"lineNum":"    6","line":""},
{"lineNum":"    7","line":"/// The Split_operator_stepper_elements class generates a constant number of"},
{"lineNum":"    8","line":"/// split-operator steps per thick element. Thin elements are assigned"},
{"lineNum":"    9","line":"/// a single step each. One or more collective effects are included per"},
{"lineNum":"   10","line":"/// step."},
{"lineNum":"   11","line":"class Split_operator_stepper_elements : public Stepper","class":"lineNoCov","hits":"0","possible_hits":"13",},
{"lineNum":"   12","line":"{"},
{"lineNum":"   13","line":"private:"},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"    int steps_per_element;"},
{"lineNum":"   16","line":"    std::shared_ptr<const CO_options> co_ops;"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"    std::vector<Step>"},
{"lineNum":"   19","line":"    apply_impl(Lattice const & lattice) const override;"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"public:"},
{"lineNum":"   22","line":""},
{"lineNum":"   23","line":"    Split_operator_stepper_elements("},
{"lineNum":"   24","line":"            CO_options const & coo = Dummy_CO_options(),"},
{"lineNum":"   25","line":"            int steps_per_element = 1 )"},
{"lineNum":"   26","line":"    : steps_per_element(steps_per_element)","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"   27","line":"    , co_ops(coo.clone())","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   28","line":"    { }","class":"lineNoCov","hits":"0","possible_hits":"10",},
{"lineNum":"   29","line":""},
{"lineNum":"   30","line":"    std::unique_ptr<Stepper> clone() const override"},
{"lineNum":"   31","line":"    { return std::make_unique<Split_operator_stepper_elements>(*this); }","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"private:"},
{"lineNum":"   34","line":""},
{"lineNum":"   35","line":"    friend class cereal::access;"},
{"lineNum":"   36","line":""},
{"lineNum":"   37","line":"    template<class Archive>"},
{"lineNum":"   38","line":"    void serialize(Archive & ar)"},
{"lineNum":"   39","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"11",},
{"lineNum":"   40","line":"        ar(cereal::base_class<Stepper>(this));","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   41","line":"        ar(steps_per_element);","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"   42","line":"        ar(co_ops);","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   43","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   44","line":"};"},
{"lineNum":"   45","line":""},
{"lineNum":"   46","line":"CEREAL_REGISTER_TYPE(Split_operator_stepper_elements)","class":"linePartCov","hits":"4","order":"695","possible_hits":"9",},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"#endif /* SPLIT_OPERATOR_STEPPER_ELEMENTS_H_ */"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 11, "covered" : 1,};
var merged_data = [];
