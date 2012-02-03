#ifndef STEP_H_
#define STEP_H_

#include <list>
#include <boost/shared_ptr.hpp>
#include "synergia/bunch/train.h"
#include "synergia/simulation/operator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/multi_diagnostics.h"
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

struct Bunch_means
{
    double x_mean;
    double y_mean;
    double z_mean;
    double realnum;
    int bucket_index;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(x_mean);
            ar & BOOST_SERIALIZATION_NVP(y_mean);
            ar & BOOST_SERIALIZATION_NVP(z_mean);
            ar & BOOST_SERIALIZATION_NVP(realnum);
            ar & BOOST_SERIALIZATION_NVP(bucket_index);
        }
};

class Step
{
private:
    Operators operators;
    std::list<double > time_fractions;
    double length;
    /// a list which contains informations about prevoius turns, the last element is the last (the earliest) turn stored
    /// it is updated at every step where impedance kick is applied
    /// the element is a vector of size num_bunches. The vector elements correspund to different bunches
    std::list< std::vector<Bunch_means> > stored_vbunches;

public:
    Step(double length);
    // Default constructor for serialization use only
    Step();
    void
    append(Operator_sptr operator_sptr, double time_fraction);
    void
    append(Operators const& operators, double time_fraction);
    virtual void
    apply(Bunch & bunch);
    virtual void
    apply(Bunch & bunch, Multi_diagnostics & diagnostics);
    virtual void
    apply(Bunch_with_diagnostics_train & bunch_diag_train);
    Operators const&
    get_operators() const;
    Operators &
    get_operators();
    std::list<double> const&
    get_time_fractions() const;
    double
    get_length() const;
    std::list< std::vector<Bunch_means> >  const& get_stored_vbunches() const;
    virtual void
    print(int index) const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(operators);
            ar & BOOST_SERIALIZATION_NVP(time_fractions);
            ar & BOOST_SERIALIZATION_NVP(length);
            ar & BOOST_SERIALIZATION_NVP(stored_vbunches);
        }
};

typedef boost::shared_ptr<Step > Step_sptr;
typedef std::list<Step_sptr > Steps;

#endif /* STEP_H_ */
