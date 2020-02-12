
#include "synergia/bunch/diagnostics_worker.h"
#include "synergia/bunch/bunch.h"

std::string Diagnostics_worker::type() const 
{ 
    return diag->type(); 
}

void Diagnostics_worker::update(Bunch const& bunch)
{ 
    diag->update(bunch); 
    diag->reduce(bunch.get_comm(), 
            diag_file.get_file().master_rank());
}

void Diagnostics_worker::write()
{ 
    diag_file.open_file();
    diag->write(diag_file.get_file()); 
    diag_file.finish_write();
}


