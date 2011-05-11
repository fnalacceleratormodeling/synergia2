#include "command_line_arg.h"

Command_line_arg::Command_line_arg(const char * char_ptr) :
    str(char_ptr)
{
}

bool
Command_line_arg::is_equal_pair() const
{
    return str.find('=') != std::string::npos;
}

std::string
Command_line_arg::get_lhs() const
{
    return str.substr(0, str.find('='));
}

std::string
Command_line_arg::get_rhs() const
{
    return str.substr(str.find('=') + 1, std::string::npos);
}
