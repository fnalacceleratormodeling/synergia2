#ifndef ERROR_HANDLER_MA_PARSE_H
#define ERROR_HANDLER_MA_PARSE_H

#include "madx.h"
#include "mx_tree.h"
#include <boost/any.hpp>

namespace synergia
{
  // parse mad8 string
  bool parse_madx( string_t const & string, MadX & mx );

  // parse mad8 file
  bool parse_madx_file( string_t const & fname, MadX & mx );
  bool parse_madx_file( string_t const & fname, mx_tree & doc );

  bool parse_madx_tree( string_t const & s, mx_tree & doc );

  bool parse_expression( std::string const & s, double & result );
}

#endif
