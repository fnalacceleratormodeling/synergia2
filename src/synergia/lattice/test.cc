#include <iostream>
#include "lattice/mx_parse.h"

int
main()
{
    synergia::MadX mx;
    synergia::parse_madx_file("ps_lattice.madx", mx);

    std::cout << "line_count = " << mx.line_count() << std::endl;
    std::cout << "sequence_count = " << mx.sequence_count() << std::endl;
    std::cout << "element count in ps sequence = "
            << mx.sequence("ps").element_count() << std::endl;
    std::cout << "ps sequence length = " << mx.sequence("ps").length()
            << std::endl;

    for (int i = 0; i < mx.sequence("ps").element_count(); ++i) {
        std::string name(mx.sequence("ps").element(i, false).label());
        if ( name == "" ) {
            name = mx.sequence("ps").element(i, true).label();
        }
        std::string type(mx.sequence("ps").element(i, true).name());
        std::cout << "element " << i << "> "
                << name << "--"
                << type;
        std::vector<string_t > attribute_names(
                mx.sequence("ps").element(i).attribute_names());
        for (std::vector<string_t >::iterator it = attribute_names.begin();
                it != attribute_names.end(); ++it){
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
    std::cout << "jfa: 886 attribute_count = "
            << mx.sequence("ps").element(886).attribute_count() << std::endl;
    std::cout << "jfa: 886 angle = "
            << mx.sequence("ps").element(886).attribute_as_number("angle")
            << std::endl;
    std::cout << "success" << std::endl;
    return 0;
}

