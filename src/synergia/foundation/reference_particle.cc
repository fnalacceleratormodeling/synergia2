#include "reference_particle.h"
#include "synergia/utils/floating_point.h"

Reference_particle::Reference_particle() :
    charge(0), four_momentum(1.0, 1.0), state(boost::extents[6]),
            repetition(0), s(0), s_n(0)
{

}

Reference_particle::Reference_particle(int charge, double mass,
        double total_energy) :
    charge(charge), four_momentum(mass, total_energy),
            state(boost::extents[6]), repetition(0), s(0), s_n(0)
{
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(int charge,
        Four_momentum const & four_momentum_in) :
    charge(charge), four_momentum(four_momentum_in), state(boost::extents[6]),
            repetition(0), s(0), s_n(0)
{
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(int charge,
        Four_momentum const & four_momentum_in, Const_MArray1d_ref state) :
    charge(charge), four_momentum(four_momentum_in), state(state),
            repetition(0), s(0), s_n(0)
{
}

Reference_particle::Reference_particle(Lsexpr const& lsexpr) :
      state(boost::extents[6])
    , repetition(0)
    , s(0)
    , s_n(0)
{
    bool found_charge(false), found_four_momentum(false);
    for(Lsexpr::const_iterator_t it = lsexpr.begin();
        it != lsexpr.end(); ++it) {
        if(it->is_labeled()) {
            if(it->get_label() == "charge") {
                charge = it->get_int();
                found_charge = true;
            } else if(it->get_label() == "four_momentum") {
                four_momentum = Four_momentum(*it);
                found_four_momentum = true;
            } else if(it->get_label() == "repetition") {
                repetition = it->get_int();
            } else if(it->get_label() == "s") {
                s = it->get_double();
            } else if(it->get_label() == "s_n") {
                s_n = it->get_double();
            } else if(it->get_label() == "state") {
                std::vector<double> state_vector(it->get_double_vector());
                if (state_vector.size() != 6) {
                    throw std::runtime_error("Reference_particle from lsexpr: state vector must have length 6");
                }
                std::copy(state_vector.begin(), state_vector.end(), state.begin());
            } else {
                throw std::runtime_error("Reference_particle: Lsexpr unknown label \"" +
                                         it->get_label() +
                                         "\"");
            }
        }
    }
    if(!found_charge) {
        throw std::runtime_error("Reference_particle: Lsexpr missing charge");
    }
    if(!found_four_momentum) {
        throw std::runtime_error("Reference_particle: Lsexpr missing four_momentum");
    }
}

Lsexpr Reference_particle::as_lsexpr() const
{
    Lsexpr retval;

    Lsexpr four_momentum_lsexpr(four_momentum.as_lsexpr());
    four_momentum_lsexpr.set_label("four_momentum");
    retval.push_back(four_momentum_lsexpr);

    retval.push_back(Lsexpr(charge, "charge"));

    if(repetition != 0) {
        retval.push_back(Lsexpr(repetition, "repetition"));
    }

    if(s != 0) {
        retval.push_back(Lsexpr(s, "s"));
    }

    if(s_n != 0) {
        retval.push_back(Lsexpr(s_n, "s_n"));
    }

    if((state[0] != 0) || (state[1] != 0) || (state[2] != 0) ||
       (state[3] != 0) || (state[4] != 0) || (state[5] != 0)) {
        std::vector<double> state_vector(6);
        std::copy(state.begin(), state.end(), state_vector.begin());
        retval.push_back(Lsexpr(state_vector, "state"));
    }

    return retval;
}

void
Reference_particle::set_four_momentum(Four_momentum const & four_momentum)
{
    this->four_momentum = four_momentum;
}

void
Reference_particle::set_state(Const_MArray1d_ref state)
{
    this->state = state;
}

void
Reference_particle::set_total_energy(double total_energy)
{
    four_momentum.set_total_energy(total_energy);
}

void
Reference_particle::increment_trajectory(double length)
{
    s_n += length;
}

void
Reference_particle::start_repetition()
{
    if (s == 0.0) {
        s = s_n;
    }
    if (s_n > 0.0) {
        repetition += 1;
    }
    s_n = 0.0;
}

void
Reference_particle::set_trajectory(int repetition, double repetition_length,
        double s)
{
    this->repetition = repetition;
    this->s = repetition_length;
    this->s_n = s;
}

int
Reference_particle::get_charge() const
{
    return charge;
}

Four_momentum const &
Reference_particle::get_four_momentum() const
{
    return four_momentum;
}

Const_MArray1d_ref
Reference_particle::get_state() const
{
    return state;
}

double
Reference_particle::get_gamma() const
{
    return four_momentum.get_gamma();
}

double
Reference_particle::get_beta() const
{
    return four_momentum.get_beta();
}

double
Reference_particle::get_momentum() const
{
    return four_momentum.get_momentum();
}

double
Reference_particle::get_total_energy() const
{
    return four_momentum.get_total_energy();
}

double
Reference_particle::get_s() const
{
    return repetition * s + s_n;
}

double
Reference_particle::get_s_n() const
{
    return s_n;
}

int
Reference_particle::get_repetition() const
{
    return repetition;
}

double
Reference_particle::get_repetition_length() const
{
    return s;
}

bool
Reference_particle::equal(Reference_particle const& reference_particle,
        double tolerance) const
{
    if (charge != reference_particle.get_charge()) {
        return false;
    }
    if (!four_momentum.equal(reference_particle.get_four_momentum(), tolerance)) {
        return false;
    }
    for (int i = 0; i < 6; ++i) {
        if (!floating_point_equal(state[i], reference_particle.get_state()[i],
                tolerance)) {
            return false;
        }
    }
    return true;
}

template<class Archive>
    void
    Reference_particle::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(charge)
                & BOOST_SERIALIZATION_NVP(four_momentum)
                & BOOST_SERIALIZATION_NVP(state)
                & BOOST_SERIALIZATION_NVP(repetition)
                & BOOST_SERIALIZATION_NVP(s)
                & BOOST_SERIALIZATION_NVP(s_n);
    }

template
void
Reference_particle::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Reference_particle::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Reference_particle::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Reference_particle::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
