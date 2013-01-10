#ifndef ERROR_HANDLER_MA_PARSE_H
#define ERROR_HANDLER_MA_PARSE_H

#include "mad8.h"

namespace synergia
{
  // parse a single mad8 statement
  bool parse_mad8_statement ( string_t const & s, Mad8 & m8 );

  // parse mad8 string
  bool parse_mad8( string_t const & string, Mad8 & m8 );

  // parse mad8 file
  bool parse_mad8_file( string_t const & fname, Mad8 & m8 );
}

#endif
