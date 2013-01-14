#include "madx_reader.h"
#include "mx_parse.h"
#include "synergia/utils/digits.h"

#include <stdexcept>
#include <algorithm>

using namespace synergia;

MadX_reader::MadX_reader()
{
}

void
MadX_reader::parse(std::string const& string)
{
    madx_sptr = boost::shared_ptr<MadX >(new MadX());
    synergia::parse_madx(string, *madx_sptr);
}

void
MadX_reader::parse_file(std::string const& filename)
{
    madx_sptr = boost::shared_ptr<MadX >(new MadX());
    synergia::parse_madx_file(filename, *madx_sptr);
}

std::vector<std::string >
MadX_reader::get_line_names()
{
    if (!madx_sptr) {
        throw std::runtime_error(
                "MadX_reader::get_line_names: nothing has been parsed");
    }
    std::vector < std::string > retval(madx_sptr->line_labels());
    std::vector < std::string > sequence_labels(madx_sptr->sequence_labels());
    for (std::vector<std::string >::const_iterator it = sequence_labels.begin();
            it != sequence_labels.end(); ++it) {
        retval.push_back(*it);
    }
    return retval;
}

Lattice
MadX_reader::get_lattice(std::string const& line_name)
{
    if (!madx_sptr) {
        throw std::runtime_error(
                "MadX_reader::get_lattice: nothing has been parsed");
    }

    Lattice lattice(line_name);

    bool found_line(false), found_sequence(false);

    std::vector < std::string > lines(madx_sptr->line_labels());
    std::vector<std::string >::const_iterator it = std::find(lines.begin(),
            lines.end(), line_name);
    if (it != lines.end()) {
        found_line = true;
    } else {
        std::vector < std::string > sequences(madx_sptr->sequence_labels());
        std::vector<std::string >::const_iterator it = std::find(
                sequences.begin(), sequences.end(), line_name);
        if (it != lines.end()) {
            found_sequence = true;
        } else {
            throw std::runtime_error(
                    "cannot find line with label " + line_name);
        }
    }
    if (found_line) {
        // to be completed
        throw std::runtime_error(
                "MadX_reader::get_lattice does not currently handle lines");
    }
    if (found_sequence) {
        MadX_sequence sequence(madx_sptr->sequence(line_name));
        double total_length = sequence.length();
        double current_pos = 0.0;
        const double min_drift_length = 1.0e-9;
        int drift_count = 0;
        int drift_digits = digits(sequence.element_count());
        for (int i = 0; i < sequence.element_count(); ++i) {
            double at = sequence.element(i).attribute_as_number("at");
            double drift_length = at - current_pos;
            if (drift_length > min_drift_length) {
                std::stringstream name_stream;
                name_stream << "auto_drift";
                name_stream << "_";
                name_stream << std::setw(drift_digits);
                name_stream << std::setfill('0');
                name_stream << drift_count;
                Lattice_element drift("drift", name_stream.str());
                drift.set_double_attribute("l", drift_length, false);
                lattice.append(drift);
            }
            std::string name(sequence.element(i, false).label());
            if (name == "") {
                name = sequence.element(i, true).label();
            }
            std::string type(sequence.element(i, true).name());
            Lattice_element element(type, name);
            std::vector<string_t > attribute_names(
                    sequence.element(i).attribute_names());
            for (std::vector<string_t >::iterator it = attribute_names.begin();
                    it != attribute_names.end(); ++it) {
                element.set_double_attribute(*it,
                        sequence.element(i).attribute_as_number(*it));
            }
            lattice.append(element);
            current_pos = at + element.get_length();
        }
    }

    return lattice;
}

Lattice
MadX_reader::get_lattice(std::string const& line_name,
        std::string const& filename)
{
    parse_file(filename);
    return get_lattice(line_name);
}

MadX_reader::~MadX_reader()
{

}
