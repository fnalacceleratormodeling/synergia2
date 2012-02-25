#include "resume.h"

Resume::Resume(std::string const& checkpoint_dir) :
    checkpoint_dir(checkpoint_dir), propagator()
{
    remove_serialization_directory();
    symlink_serialization_directory(checkpoint_dir);
    binary_load(
            propagator,
            get_combined_path(checkpoint_dir,
                    Propagator::propagator_archive_name).c_str());
    unlink_serialization_directory();
}

void
Resume::set_checkpoint_period(int period)
{
    propagator.set_checkpoint_period(period);
}

int
Resume::get_checkpoint_period() const
{
    return propagator.get_checkpoint_period();
}

void
Resume::set_new_checkpoint_dir(std::string const& directory_name)
{
    propagator.set_checkpoint_dir(directory_name);
}

std::string const&
Resume::get_new_checkpoint_dir() const
{
    return propagator.get_checkpoint_dir();
}

void
Resume::propagate()
{
    propagator.resume(checkpoint_dir);
}

void
Resume::propagate(int max_turns)
{
    propagator.resume(checkpoint_dir, max_turns);
}

Resume::~Resume()
{

}
