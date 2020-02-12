
#include "synergia/bunch/diagnostics_worker.h"
#include "synergia/bunch/bunch.h"

std::string Diagnostics_worker::type() const 
{ 
    return diag->type(); 
}

void Diagnostics_worker::update(Bunch const& bunch)
{ 
    diag->update(bunch); 
    diag->reduce(bunch.get_comm(), writer.get_file().master_rank());
}

void Diagnostics_worker::write()
{ 
    writer.open_file();
    diag->write(writer.get_file()); 
    writer.finish_write();
}


