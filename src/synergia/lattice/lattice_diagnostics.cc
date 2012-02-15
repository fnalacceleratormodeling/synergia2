#include "lattice_diagnostics.h"
#include "synergia/utils/hdf5_writer.h"

Lattice_diagnostics::Lattice_diagnostics(Lattice_sptr lattice_sptr,
        std::string const& filename, std::string const& attribute) :
    Generalized_diagnostics("lattice_diagnostics"), lattice_sptr(lattice_sptr),
            filename(filename), attribute(attribute), default_value(0),
            reduce(true), reduce_op(MPI_SUM), initial_lattice_size(0),
            write_helper(filename, true, Commxx(MPI_COMM_WORLD)), writer(0),
            first_time(true)

{
}

void
Lattice_diagnostics::set_default_value(double value)
{
    default_value = value;
}

double
Lattice_diagnostics::get_default_value() const
{
    return default_value;
}

void
Lattice_diagnostics::set_reduce(bool reduce)
{
    this->reduce = reduce;
}

bool
Lattice_diagnostics::get_reduce() const
{
    return reduce;
}

void
Lattice_diagnostics::set_reduce_op(MPI_Op op)
{
    reduce_op = op;
}

MPI_Op
Lattice_diagnostics::get_reduce_op() const
{
    return reduce_op;
}

bool
Lattice_diagnostics::is_serial() const
{
    return true;
}

void
Lattice_diagnostics::update()
{
}

void
Lattice_diagnostics::write()
{
    int lattice_size = lattice_sptr->get_elements().size();
    if (first_time) {
        initial_lattice_size = lattice_size;
    }
    if (lattice_size != initial_lattice_size) {
        throw std::runtime_error(
                "Lattice_diagnostics::write: Number of lattice elements has changed since the first write");
    }
    if (first_time) {
        if (write_helper.write_locally()) {
            MArray1d ss(boost::extents[lattice_size]);
            int index = 0;
            double s = 0.0;
            for (Lattice_elements::const_iterator it =
                    lattice_sptr->get_elements().begin(); it
                    != lattice_sptr->get_elements().end(); ++it) {
                s += (*it)->get_length();
                ss[index] = s;
                ++index;
            }
            write_helper.get_hdf5_file()->write(ss, "s");
            writer = new Hdf5_serial_writer<MArray1d_ref > (
                    write_helper.get_hdf5_file(), attribute);
        }
        first_time = false;
    }
    MArray1d values(boost::extents[lattice_size]);
    int index = 0;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        double value = default_value;
        if ((*it)->has_double_attribute(attribute)) {
            value = (*it)->get_double_attribute(attribute);
        }
        values[index] = value;
        ++index;
    }
    MArray1d reduced_values(boost::extents[lattice_size]);
    if (reduce) {
        MPI_Reduce((void*) values.origin(), (void *) reduced_values.origin(),
                values.num_elements(), MPI_DOUBLE, reduce_op,
                write_helper.get_writer_rank(), MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < lattice_size; ++i) {
            reduced_values[i] = values[i];
        }
    }
    if (write_helper.write_locally()) {
        writer->append(reduced_values);
        write_helper.finish_write();
    }
}

Lattice_diagnostics::~Lattice_diagnostics()
{
    if (!first_time) {
        if (write_helper.write_locally()) {
            delete writer;
        }
    }
}
