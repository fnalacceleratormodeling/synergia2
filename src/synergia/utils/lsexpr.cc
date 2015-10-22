#include "lsexpr.h"

#include <iostream>
Lsexpr read_lsexpr(std::istream & stream)
{
    Lsexpr retval;
    bool in_str(false);
    Lsexpr * current = &retval;
    Lsexpr * parent = &retval;
    char c(stream.get());
    while(!stream.eof()) {
        if(c == '{' && !in_str) {
            Lsexpr empty;
            current->sequence.push_back(empty);
            parent = current;
            current = &(current->sequence.back());
        } else if(c == '}' && !in_str) {
            if(!current->atom.empty()) {
                current->is_atom = true;
            }
            Lsexpr empty;
            current = &(parent->sequence.back());
            current->sequence.push_back(empty);
        } else if((c == ',') && !in_str) {
            if(!current->atom.empty()) {
                current->is_atom = true;
            }
            Lsexpr empty;
            current = parent;
            current->sequence.push_back(empty);
            current = &(current->sequence.back());
        } else if(c=='"') {
            in_str = !in_str;
        } else if((c == ' ') || (c == '\n') || (c == '\t')) {
            // pass
        } else if((c == ':' && !in_str)) {
            current->has_label = true;
            current->label = current->atom;
            current->atom = "";
        } else {
            current->atom += c;
        }
        c = stream.get();
    }
    if(!current->atom.empty()) {
        current->is_atom = true;
    }
    return retval;
}

namespace {
void write_lsexpr_internal(Lsexpr const& lsexpr, std::ostream & stream,
                  int indent, bool first)
{
    if(lsexpr.has_label) {
        stream << lsexpr.label << ": ";
    }
    if(lsexpr.is_atom) {
        stream << lsexpr.atom;
    } else {
        if(!first || lsexpr.has_label) {
            if(!lsexpr.is_atom) {
                indent += 4;
            }
            stream << "\n" << std::string(indent, ' ');
        }
        stream << "{";
        for(std::vector<Lsexpr>::const_iterator it = lsexpr.sequence.begin();
            it != lsexpr.sequence.end(); ++it) {
            write_lsexpr_internal(*it, stream, indent, false);
            if ((it+1) != lsexpr.sequence.end()) {
                stream << ", ";
            }
        }
        stream << "}";
        if(!first || lsexpr.has_label) {
            indent -= 4;
            stream << "\n";
            if (indent > 0) {
                stream << std::string(indent, ' ');
            }
        }
    }
}
}

void write_lsexpr(Lsexpr const& lsexpr, std::ostream & stream)
{
    write_lsexpr_internal(lsexpr, stream, 0, true);
}
