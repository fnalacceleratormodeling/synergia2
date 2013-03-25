#ifndef ERROR_HANDLER_MA_PARSE_H
#define ERROR_HANDLER_MA_PARSE_H

#include "madx.h"
#include "mx_tree.h"
#include <boost/any.hpp>

namespace synergia
{
  // parse into a MadX object
  bool parse_madx( string_t const & string, MadX & mx );
  bool parse_madx_file( string_t const & fname, MadX & mx );

  // parse into an intermediate mx_tree object
  bool parse_int_madx( string_t const & s, mx_tree & doc );
  bool parse_int_madx_file( string_t const & fname, mx_tree & doc );

  // parse a const math expression
  bool parse_expression( std::string const & s, double & result );
}

#endif
