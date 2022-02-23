#include "synergia/lattice/lattice_tree.h"
#include "synergia/lattice/mx_parse.h"

void
Lattice_tree::set_variable(std::string const& name, double val)
{
    mx.insert_variable(name, synergia::mx_expr(val));
}

void
Lattice_tree::set_variable(std::string const& name, std::string const& val)
{
    synergia::mx_expr expr;
    synergia::parse_expression(val, expr);

    mx.insert_variable(name, expr);
}

void
Lattice_tree::print() const
{

}


