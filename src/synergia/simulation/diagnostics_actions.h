#ifndef DIAGNOSTICS_ACTIONS_H_
#define DIAGNOSTICS_ACTIONS_H_

#include <list>

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/simulation/stepper.h"

class Diagnostics_actions
{
public:
    struct Periodic
    {
        int period;
        Diagnostics_sptr diagnostics_sptr;
        Periodic(int period, Diagnostics_sptr diagnostics_sptr) :
            period(period), diagnostics_sptr(diagnostics_sptr)
        {
        }
        Periodic()
        {
        }
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version);
    };
    typedef std::list<Periodic > Periodics;

    typedef std::list<int > Numbers;

    struct Listed
    {
        Numbers numbers;
        Diagnostics_sptr diagnostics_sptr;
        Listed(Numbers const& numbers, Diagnostics_sptr diagnostics_sptr) :
            numbers(numbers), diagnostics_sptr(diagnostics_sptr)
        {
        }
        Listed()
        {
        }
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version);
    };
    typedef std::list<Listed > Listeds;

    struct Periodic_listed
    {
        int turn_period;
        Numbers step_numbers;
        Diagnostics_sptr diagnostics_sptr;
        Periodic_listed(int turn_period, Numbers const& step_numbers,
                Diagnostics_sptr diagnostics_sptr) :
            turn_period(turn_period), step_numbers(step_numbers),
                    diagnostics_sptr(diagnostics_sptr)
        {
        }
        Periodic_listed()
        {
        }
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version);
    };
    typedef std::list<Periodic_listed > Periodic_listeds;

private:
    Bunch_sptr bunch_sptr;
    bool have_bunch_sptr_;
    Periodics per_turn_periodic, per_step_periodic;
    Listeds per_turn_listed;
    Periodic_listeds per_step_periodic_listed;

    void
    update_and_write_periodics(Periodics & periodics, int num);
    void
    update_and_write_listeds(Listeds & listeds, int num);
    void
    update_and_write_periodic_listeds(Periodic_listeds & periodic_listeds,
            int step, int turn);

public:
    Diagnostics_actions();
    virtual void
    set_bunch_sptr(Bunch_sptr bunch_sptr);
    virtual bool
    have_bunch_sptr() const;
    virtual Bunch_sptr
    get_bunch_sptr();
    virtual void
    add_per_turn(Diagnostics_sptr diagnostics_sptr, int period = 1);
    virtual void
    add_per_turn(Diagnostics_sptr diagnostics_sptr,
            std::list<int > const& turn_numbers);
    virtual void
    add_per_step(Diagnostics_sptr diagnostics_sptr, int period = 1);
    virtual void
    add_per_step(Diagnostics_sptr diagnostics_sptr,
            std::list<int > const& step_numbers, int turn_period = 1);
    virtual void
    first_action(Stepper & stepper, Bunch & bunch);
    virtual void
    turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num);
    virtual void
    step_end_action(Stepper & stepper, Step & step, Bunch & bunch,
            int turn_num, int step_num);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Diagnostics_actions();
};

typedef boost::shared_ptr<Diagnostics_actions >
        Diagnostics_actions_sptr;
#endif /* DIAGNOSTICS_ACTIONS_H_ */
