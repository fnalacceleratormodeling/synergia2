var data = {lines:[
{"lineNum":"    1","line":""},
{"lineNum":"    2","line":"#include \"synergia/bunch/diagnostics_worker.h\""},
{"lineNum":"    3","line":"#include \"synergia/bunch/bunch.h\""},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":"std::string Diagnostics_worker::type() const"},
{"lineNum":"    6","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"    7","line":"    return diag->type();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"    8","line":"}"},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"void Diagnostics_worker::update(Bunch const& bunch)"},
{"lineNum":"   11","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   12","line":"    diag->update(bunch);"},
{"lineNum":"   13","line":"    diag->reduce(bunch.get_comm(),"},
{"lineNum":"   14","line":"            diag_file.get_file().master_rank());","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   15","line":"}"},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"void Diagnostics_worker::write()"},
{"lineNum":"   18","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   19","line":"    diag_file.open_file();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   20","line":"    diag->write(diag_file.get_file());","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   21","line":"    diag_file.finish_write();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   22","line":"}"},
{"lineNum":"   23","line":""},
{"lineNum":"   24","line":""},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 8, "covered" : 0,};
var merged_data = [];
