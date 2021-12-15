var data = {lines:[
{"lineNum":"    1","line":"#ifndef STRING_UTILS_H_"},
{"lineNum":"    2","line":"#define STRING_UTILS_H_"},
{"lineNum":"    3","line":"#include <string>"},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":"inline bool"},
{"lineNum":"    6","line":"false_string(std::string const& val)"},
{"lineNum":"    7","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"    8","line":"    return (val == \"false\") || (val == \"False\") || (val == \"FALSE\");","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"    9","line":"}"},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"#endif /* STRING_UTILS_H_ */"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 2, "covered" : 0,};
var merged_data = [];
