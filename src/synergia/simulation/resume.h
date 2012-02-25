#ifndef RESUME_H_
#define RESUME_H_

#include <string>
#include "synergia/simulation/propagator.h"

class Resume
{
private:
    std::string checkpoint_dir;
    Propagator propagator;
public:
    Resume(std::string const& checkpoint_dir = Propagator::default_checkpoint_dir);

    void
    set_checkpoint_period(int period);

    int
    get_checkpoint_period() const;

    void
    set_new_checkpoint_dir(std::string const& directory_name);

    std::string const&
    get_new_checkpoint_dir() const;

    void
    propagate();

    void
    propagate(int max_turns);

    ~Resume();
};

#endif /* RESUME_H_ */
