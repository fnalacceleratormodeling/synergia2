
#include "synergia/bunch/diagnostics_worker.h"
#include "synergia/bunch/bunch.h"

std::string Diagnostics_worker::type() const 
{ 
    return diag->type(); 
}

void Diagnostics_worker::update(Bunch const& bunch)
{ 
    diag->update(bunch); 
    diag->collect(bunch.get_comm(), writer.writer_rank());
}

void Diagnostics_worker::write()
{ 
    if (writer.write_locally())
    { 
        writer.open_file();
        diag->write(writer.get_file()); 
        writer.finish_write(diag->serial());
    }
}


