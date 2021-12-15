var data = {lines:[
{"lineNum":"    1","line":""},
{"lineNum":"    2","line":"#include \"lattice_element_processor.h\""},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":""},
{"lineNum":"    6","line":"Lattice_element"},
{"lineNum":"    7","line":"Lattice_element_processor::process(Lattice_element const & element)"},
{"lineNum":"    8","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"    9","line":"    auto t = element.get_type();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   10","line":"    Lattice_element e(element);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   11","line":""},
{"lineNum":"   12","line":"    switch(element.get_type())","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   13","line":"    {"},
{"lineNum":"   14","line":"    case element_type::drift:      drift(e); break;"},
{"lineNum":"   15","line":"    case element_type::sbend:      sbend(e); break;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   16","line":"    case element_type::quadrupole: quadrupole(e); break;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   17","line":"    case element_type::multipole:  multipole(e); break;"},
{"lineNum":"   18","line":"    case element_type::rfcavity:   rfcavity(e); break;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"    case element_type::hkicker:    hkicker(e); break;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   21","line":"    case element_type::vkicker:    vkicker(e); break;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   22","line":"    case element_type::kicker:     kicker(e); break;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   23","line":"    default: break;"},
{"lineNum":"   24","line":"    }"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"    return e;"},
{"lineNum":"   27","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   28","line":""},
{"lineNum":"   29","line":"void Lattice_element_processor::drift(Lattice_element & e)"},
{"lineNum":"   30","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   31","line":"    e.set_default_double_attribute(\"l\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   32","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   33","line":""},
{"lineNum":"   34","line":"void Lattice_element_processor::sbend(Lattice_element & e)"},
{"lineNum":"   35","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   36","line":"    e.set_default_string_attribute(\"propagator_type\", \"yoshida\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   37","line":"    e.set_default_double_attribute(\"l\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   38","line":"    e.set_default_double_attribute(\"angle\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   39","line":"    e.set_default_double_attribute(\"tilt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   40","line":"    e.set_default_double_attribute(\"k1\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   41","line":"    e.set_default_double_attribute(\"e1\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   42","line":"    e.set_default_double_attribute(\"e2\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   43","line":"    e.set_default_double_attribute(\"fint\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   44","line":"    e.set_default_double_attribute(\"fintx\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   45","line":"    e.set_default_double_attribute(\"hgap\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   46","line":"    e.set_default_double_attribute(\"k2\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   47","line":"    e.set_default_double_attribute(\"h1\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   48","line":"    e.set_default_double_attribute(\"h2\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   49","line":"    e.set_default_double_attribute(\"kicks\", 40.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   50","line":""},
{"lineNum":"   51","line":"    // possible higher order multipole components"},
{"lineNum":"   52","line":"    e.set_default_double_attribute(\"kl\", 0.0); // base strength/B-rho","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   53","line":"    e.set_default_double_attribute(\"a1\", 0.0); // skew quad","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   54","line":"    e.set_default_double_attribute(\"a2\", 0.0); // skew sextupole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   55","line":"    e.set_default_double_attribute(\"a3\", 0.0); // skew octupole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   56","line":"    e.set_default_double_attribute(\"a4\", 0.0); // skew decapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   57","line":"    e.set_default_double_attribute(\"a5\", 0.0); // skew dodecapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   58","line":"    e.set_default_double_attribute(\"a6\", 0.0); // skew tetradecapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   59","line":"    e.set_default_double_attribute(\"a7\", 0.0); // skew hexdecapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   60","line":"    e.set_default_double_attribute(\"b1\", 0.0); // quad","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   61","line":"    e.set_default_double_attribute(\"b2\", 0.0); // sextupole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   62","line":"    e.set_default_double_attribute(\"b3\", 0.0); // octopole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   63","line":"    e.set_default_double_attribute(\"b4\", 0.0); // decapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   64","line":"    e.set_default_double_attribute(\"b5\", 0.0); // dodecapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   65","line":"    e.set_default_double_attribute(\"b6\", 0.0); // tetradecapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   66","line":"    e.set_default_double_attribute(\"b7\", 0.0); // hexdecapole","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   67","line":"    e.set_default_double_attribute(\"entry_edge_kick\", 1.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   68","line":"    e.set_default_double_attribute(\"exit_edge_kick\", 1.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   69","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"void Lattice_element_processor::quadrupole(Lattice_element & e)"},
{"lineNum":"   72","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   73","line":"    e.set_default_string_attribute(\"propagator_type\", \"yoshida\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   74","line":"    e.set_default_double_attribute(\"yoshida_order\", 2.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   75","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   76","line":""},
{"lineNum":"   77","line":"void Lattice_element_processor::multipole(Lattice_element & e)"},
{"lineNum":"   78","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   79","line":"    e.set_default_double_attribute(\"tilt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   80","line":"    //e.set_default_vector_attribute(\"knl\", std::vector<double >(0));"},
{"lineNum":"   81","line":"    //e.set_default_vector_attribute(\"ksl\", std::vector<double >(0));"},
{"lineNum":"   82","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   83","line":""},
{"lineNum":"   84","line":"void Lattice_element_processor::hkicker(Lattice_element& e)"},
{"lineNum":"   85","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   86","line":"    double hk = 0.0;"},
{"lineNum":"   87","line":""},
{"lineNum":"   88","line":"    if (e.has_double_attribute(\"kick\"))","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   89","line":"        hk = e.get_double_attribute(\"kick\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   90","line":"    else if (e.has_double_attribute(\"hkick\"))","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   91","line":"        hk = e.get_double_attribute(\"hkick\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   92","line":""},
{"lineNum":"   93","line":"    // defaults"},
{"lineNum":"   94","line":"    e.set_default_double_attribute(\"l\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   95","line":"    e.set_default_double_attribute(\"tilt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   96","line":""},
{"lineNum":"   97","line":"    // forced"},
{"lineNum":"   98","line":"    e.set_double_attribute(\"hkick\", hk);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   99","line":"    e.set_double_attribute(\"vkick\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  100","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  101","line":""},
{"lineNum":"  102","line":"void Lattice_element_processor::vkicker(Lattice_element& e)"},
{"lineNum":"  103","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  104","line":"    double vk = 0.0;"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"    if (e.has_double_attribute(\"kick\"))","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  107","line":"        vk = e.get_double_attribute(\"kick\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  108","line":"    else if (e.has_double_attribute(\"vkick\"))","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  109","line":"        vk = e.get_double_attribute(\"vkick\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  110","line":""},
{"lineNum":"  111","line":"    // defaults"},
{"lineNum":"  112","line":"    e.set_default_double_attribute(\"l\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  113","line":"    e.set_default_double_attribute(\"tilt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  114","line":""},
{"lineNum":"  115","line":"    // forced"},
{"lineNum":"  116","line":"    e.set_double_attribute(\"hkick\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  117","line":"    e.set_double_attribute(\"vkick\", vk);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  118","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  119","line":""},
{"lineNum":"  120","line":"void Lattice_element_processor::kicker(Lattice_element& e)"},
{"lineNum":"  121","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  122","line":"    // default"},
{"lineNum":"  123","line":"    e.set_default_double_attribute(\"l\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  124","line":"    e.set_default_double_attribute(\"hkick\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  125","line":"    e.set_default_double_attribute(\"vkick\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  126","line":"    e.set_default_double_attribute(\"tilt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  127","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  128","line":""},
{"lineNum":"  129","line":"void Lattice_element_processor::rfcavity(Lattice_element & e)"},
{"lineNum":"  130","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  131","line":"    e.set_default_double_attribute(\"l\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  132","line":"    e.set_default_double_attribute(\"volt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  133","line":"    e.set_default_double_attribute(\"lag\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  134","line":"    e.set_default_double_attribute(\"harmon\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  135","line":"    e.set_default_double_attribute(\"shunt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  136","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  137","line":""},
{"lineNum":"  138","line":""},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":""},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 87, "covered" : 0,};
var merged_data = [];
