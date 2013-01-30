#include <iostream>
#include <stdexcept>
#include "mx_parse.h"

typedef std::vector<std::string > Strings;

void
show_beam(synergia::MadX const& mx)
{
    const int not_found = -1;
    int beam_index = not_found;
    Strings command_names(mx.commands());
    int current_index = 0;
    for (Strings::const_iterator it = command_names.begin();
            it != command_names.end(); ++it) {
        if (*it == "beam") {
            beam_index = current_index;
        }
        ++current_index;
    }
    if (beam_index > not_found) {
        synergia::MadX_command command(mx.command(beam_index));
        Strings attributes(command.attribute_names());
        std::cout << "beam attributes:";
        for (Strings::const_iterator it = attributes.begin();
                it != attributes.end(); ++it) {
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
}

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

    show_beam(mx);
    std::cout << "command_count = " << mx.command_count() << std::endl;
    if (mx.command_count() > 0) {
        Strings command_names(mx.commands());
        std::cout << "unlabeled command names:";
        for (Strings::const_iterator it = command_names.begin();
                it != command_names.end(); ++it) {
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
    std::cout << "label_count = " << mx.label_count() << std::endl;
    if (mx.label_count() > 0) {
        Strings command_labels(mx.command_labels());
        std::cout << "labeled command names:";
        for (Strings::const_iterator it = command_labels.begin();
                it != command_labels.end(); ++it) {
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
    std::cout << "command_labels = " << mx.command_count() << std::endl;
    if (mx.command_count() > 0) {
        Strings command_names(mx.command_labels());
        std::cout << "line names:";
        for (Strings::const_iterator it = command_names.begin();
                it != command_names.end(); ++it) {
            std::cout << " " << *it;
        }
        std::cout << std::endl;
    }
    std::cout << "line_count = " << mx.line_count() << std::endl;
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

