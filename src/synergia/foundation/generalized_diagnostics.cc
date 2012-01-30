#include "generalized_diagnostics.h"

Generalized_diagnostics::Generalized_diagnostics(std::string const& name) :
    name(name)
{
}

std::string const &
Generalized_diagnostics::get_name() const
{
    return name;
}

void
Generalized_diagnostics::update_and_write()
{
    update();
    write();
}
