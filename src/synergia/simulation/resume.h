#ifndef RESUME_H_
#define RESUME_H_

#include <string>
#include "synergia/simulation/propagator.h"

class Resume
{
public:
    struct Content
    {
        Bunch_sptr bunch_sptr;
        Stepper_sptr stepper_sptr;
        Lattice_sptr lattice_sptr;
        Content(Bunch_simulator * bunch_simulator_ptr,
                Stepper_sptr stepper_sptr);
    };

private:
    std::string checkpoint_dir;
    Propagator propagator;
public:
    Resume(
            std::string const& checkpoint_dir =
                    Propagator::default_checkpoint_dir);

    void
    set_checkpoint_period(int period);

    int
    get_checkpoint_period() const;

    void
    set_new_checkpoint_dir(std::string const& directory_name);

    std::string const&
    get_new_checkpoint_dir() const;

    void
    set_checkpoint_with_xml(bool with_xml);

    bool
    get_checkpoint_with_xml() const;

    void
    set_final_checkpoint(bool final_checkpoint);

    bool
    get_final_checkpoint() const;

    void
    set_concurrent_io(int max);

    int
    get_concurrent_io() const;

    /*void
    set_num_turns(int num_turns);

    int
    get_num_turns() const;
*/
    Content
    get_content();

    void
    propagate(bool new_max_turns, int max_turns, bool new_verbosity,
            int verbosity);

    ~Resume();
};

#endif /* RESUME_H_ */
