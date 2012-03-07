#ifndef STANDARD_DIAGNOSTICS_ACTIONS_H_
#define STANDARD_DIAGNOSTICS_ACTIONS_H_

#include "synergia/simulation/diagnostic_actions.h"
#include "synergia/bunch/diagnostics.h"
#include <list>

class Standard_diagnostics_actions : public Diagnostic_actions
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
            serialize(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(period);
                ar & BOOST_SERIALIZATION_NVP(diagnostics_sptr);
            }
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
            serialize(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(numbers);
                ar & BOOST_SERIALIZATION_NVP(diagnostics_sptr);
            }
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
            serialize(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(turn_period);
                ar & BOOST_SERIALIZATION_NVP(step_numbers);
                ar & BOOST_SERIALIZATION_NVP(diagnostics_sptr);
            }
    };
    typedef std::list<Periodic_listed > Periodic_listeds;

private:
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
    Standard_diagnostics_actions();
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
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostic_actions);
            ar & BOOST_SERIALIZATION_NVP(per_turn_periodic);
            ar & BOOST_SERIALIZATION_NVP(per_step_periodic);
            ar & BOOST_SERIALIZATION_NVP(per_turn_listed);
            ar & BOOST_SERIALIZATION_NVP(per_step_periodic_listed);
        }
    virtual
    ~Standard_diagnostics_actions();
};

BOOST_CLASS_EXPORT_KEY(Standard_diagnostics_actions)
;
typedef boost::shared_ptr<Standard_diagnostics_actions >
        Standard_diagnostics_actions_sptr;
#endif /* STANDARD_DIAGNOSTICS_ACTIONS_H_ */
