#include <iostream>
#include <stdexcept>
#include "mx_parse.h"

int
main(int argc, char ** argv)
{
    if (argc < 2) {
        throw std::runtime_error(
                std::string(argv[0]) + " usage: <file> [line name]");
    }
    char *filename = argv[1];
    synergia::MadX mx;
    synergia::parse_madx_file(filename, mx);

    std::cout << "line_count = " << mx.line_count() << std::endl;
    typedef std::vector<std::string > Strings;
    if (mx.line_count() > 0) {
        Strings line_names(mx.line_labels());
        std::cout << "line names:";
        for (Strings::const_iterator it = line_names.begin();
                it != line_names.end(); ++it) {
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
    std::cout << "sequence_count = " << mx.sequence_count() << std::endl;
    if (mx.sequence_count() > 0) {
        Strings sequence_names(mx.sequence_labels());
        std::cout << "sequence names:";
        for (Strings::const_iterator it = sequence_names.begin();
                it != sequence_names.end(); ++it) {
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
    if (argc > 2) {
        std::string linename(argv[2]);
        std::cout << "element count in '" << linename << "' sequence = "
        << mx.sequence(linename).element_count() << std::endl;
        std::cout << "fodo sequence length = " << mx.sequence(linename).length()
        << std::endl;

        for (size_t i = 0; i < mx.sequence(linename).element_count(); ++i) {
            std::string name(mx.sequence(linename).element(i, false).label());
            if (name == "") {
                name = mx.sequence(linename).element(i, true).label();
            }
            std::string type(mx.sequence(linename).element(i, true).name());
            std::cout << "element " << i << "> " << name << "--" << type;
            std::vector<string_t > attribute_names(
                    mx.sequence(linename).element(i).attribute_names());
            for (std::vector<string_t >::iterator it = attribute_names.begin();
                    it != attribute_names.end(); ++it) {
                std::cout << " " << *it;
            }
            std::cout << std::endl;
        }
//        std::cout << "jfa: 886 attribute_count = "
//                << mx.sequence(linename).element(886).attribute_count() << std::endl;
//        std::cout << "jfa: 886 angle = "
//                << mx.sequence(linename).element(886).attribute_as_number("angle")
//                << std::endl;
    }
    std::cout << "success" << std::endl;
    return 0;
}

