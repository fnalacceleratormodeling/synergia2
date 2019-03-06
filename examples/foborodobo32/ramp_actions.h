#ifndef __RAMP_ACTIONS__
#define __RAMP_ACTIONS__
class Ramp_actions: public Propagate_actions {
public:
    Ramp_actions(); // for serialization use
    Ramp_actions(int initdata);
    void turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num);
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version);
    ~Ramp_actions();
private:
    int data;
};

BOOST_CLASS_EXPORT_KEY(Ramp_actions)

#endif // __RAMP_ACTIONS__
