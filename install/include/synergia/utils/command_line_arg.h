#ifndef COMMAND_LINE_ARG_H_
#define COMMAND_LINE_ARG_H_
#include <string>
#include <sstream>
#include <stdexcept>

class Command_line_arg
{
private:
    std::string str;
public:
    Command_line_arg(const char * char_ptr);
    bool
    is_equal_pair() const;
    std::string
    get_lhs() const;
    std::string
    get_rhs() const;
    template<typename T>
        T
        extract_value() const
        {
            std::stringstream str_stream(get_rhs());
            T retval;
            str_stream >> retval;
            return retval;
        }
};

template<>
    bool
    Command_line_arg::extract_value<bool >() const
    {
        bool retval;
        if ((get_rhs() == "True") || (get_rhs() == "true")
                || (get_rhs() == "t") || (get_rhs() == "1")) {
            retval = true;
        } else if ((get_rhs() == "False") || (get_rhs() == "false")
                || (get_rhs() == "f") || (get_rhs() == "0") || (get_rhs()
                == "nil")) {
            retval = false;
        } else {
            throw std::runtime_error("Cannot convert " + get_rhs() + " to bool");
        }
        return retval;
    }

#endif /* COMMAND_LINE_ARG_ */
