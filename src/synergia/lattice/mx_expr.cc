#include "mx_expr.h"
#include "madx.h"

#include <stdexcept>

using namespace synergia;

double 
  mx_calculator::operator()(double val) const
{
  return val;
}

double 
  mx_calculator::operator()(std::string const & ref) const
{
  if( mx==NULL )
    if( has_def )  return def;
    else throw std::runtime_error("Unable to locate reference " + ref);

  return has_def ? mx->variable_as_number(ref, def) 
                 : mx->variable_as_number(ref);
}

double 
  mx_calculator::operator()(string_pair_t const & ref) const
{
  if( mx==NULL )
    if( has_def ) return def;
    else throw std::runtime_error("Unable to locate reference " + ref.first + "->" + ref.second);

  return has_def ? mx->command(ref.first).attribute_as_number(ref.second, def)
                 : mx->command(ref.first).attribute_as_number(ref.second);
}

double 
  mx_calculator::operator()(nop_t const & n) const
{
  return boost::apply_visitor(*this, n.expr);
}

double 
  mx_calculator::operator()(uop_t const & u) const
{
  return u.func( boost::apply_visitor(*this, u.param) );
}

double 
  mx_calculator::operator()(bop_t const & b) const
{
  return b.func( boost::apply_visitor(*this, b.lhs)
               , boost::apply_visitor(*this, b.rhs) );
}


