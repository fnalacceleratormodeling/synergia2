#include <stdexcept>
#include <sstream>
#include "mad8_adaptors.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <beamline/beamline_elements.h>
#include <beamline/YoshidaPropagator.h>
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

Marker_mad8_adaptor::Marker_mad8_adaptor()
{
}

Chef_elements
Marker_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(new marker(lattice_element.get_name().c_str()));
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Marker_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Marker_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Marker_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Marker_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Marker_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Marker_mad8_adaptor::~Marker_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Marker_mad8_adaptor)

Drift_mad8_adaptor::Drift_mad8_adaptor()
{
}

Chef_elements
Drift_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(
            new drift(lattice_element.get_name().c_str(),
                    lattice_element.get_length()));
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Drift_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Drift_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Drift_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Drift_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Drift_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Drift_mad8_adaptor::~Drift_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Drift_mad8_adaptor)

Sbend_mad8_adaptor::Sbend_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("angle", 0.0);
    get_default_element().set_double_attribute("k1", 0.0);
    get_default_element().set_double_attribute("e1", 0.0);
    get_default_element().set_double_attribute("e2", 0.0);
    get_default_element().set_double_attribute("k2", 0.0);
    get_default_element().set_double_attribute("h1", 0.0);
    get_default_element().set_double_attribute("h2", 0.0);
    get_default_element().set_double_attribute("hgap", 0.0);
    get_default_element().set_double_attribute("fint", 0.0);
    get_default_element().set_double_attribute("k3", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
    get_default_element().set_double_attribute("kicks", 40.0);
    // possible higher order multipole components
    get_default_element().set_double_attribute("kl", 0.0); // base strength/B-rho
    get_default_element().set_double_attribute("a1", 0.0); // skew quad
    get_default_element().set_double_attribute("a2", 0.0); // skew sextupole
    get_default_element().set_double_attribute("a3", 0.0); // skew octupole
    get_default_element().set_double_attribute("a4", 0.0); // skew decapole
    get_default_element().set_double_attribute("a5", 0.0); // skew dodecapole
    get_default_element().set_double_attribute("a6", 0.0); // skew tetradecapole
    get_default_element().set_double_attribute("a7", 0.0); // skew hexdecapole
    get_default_element().set_double_attribute("b1", 0.0); // quad
    get_default_element().set_double_attribute("b2", 0.0); // sextupole
    get_default_element().set_double_attribute("b3", 0.0); // octopole
    get_default_element().set_double_attribute("b4", 0.0); // decapole
    get_default_element().set_double_attribute("b5", 0.0); // dodecapole
    get_default_element().set_double_attribute("b6", 0.0); // tetradecapole
    get_default_element().set_double_attribute("b7", 0.0); // hexdecapole
}

Chef_elements
Sbend_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");
    double k1 = lattice_element.get_double_attribute("k1");
    double k2 = lattice_element.get_double_attribute("k2");
    double k3 = lattice_element.get_double_attribute("k3");
    double tilt = lattice_element.get_double_attribute("tilt");
    int kicks = lattice_element.has_double_attribute("kicks") ?
            floor(lattice_element.get_double_attribute("kicks")) : 40;

    bool simple = ((k1 == 0.0) && (k2 == 0.0) && (k3 == 0.0));

    alignmentData aligner;
    aligner.xOffset = 0.0;
    aligner.yOffset = 0.0;
    aligner.tilt = tilt;

    double ak[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently a0-a7
    double bk[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently b0-b7

    // thinpole strengths
    string a_attr_list[] = { "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7" };
    string b_attr_list[] = { "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7" };

    bool has_multipoles = false;
    int highest_order = 0;
    // find any possible multipole moments
    for (int moment = 0; moment < 8; ++moment) {
        if (lattice_element.has_double_attribute(a_attr_list[moment])) {
            ak[moment] = lattice_element.get_double_attribute(
                    a_attr_list[moment]);
            // no point in setting using multipole machinery ff the attributes are 0
            if (ak[moment] != 0.0) {
                has_multipoles = true;
                highest_order = moment;
            }
        }
        if (lattice_element.has_double_attribute(b_attr_list[moment])) {
            bk[moment] = lattice_element.get_double_attribute(
                    b_attr_list[moment]);
            if (bk[moment] != 0.0) {
                has_multipoles = true;
                highest_order = moment;
            }
        }
    }

    // being not simple or having a tilt precludes using multipoles
    if ((!simple || (tilt != 0.0)) && has_multipoles) {
        throw runtime_error(
                "shouldn't use k1,k2,tilt and multipoles in one sbend");
    }

    if (simple) {

        ElmPtr elm(
                new sbend(lattice_element.get_name().c_str(), length,
                        brho * angle / length, angle, e1, e2));
        elm->setTag("SBEND");
        if (tilt != 0.0) elm->setAlignment(aligner);

        // if there are no multipoles, I'm done.
        if (!has_multipoles) {
            retval.push_back(elm);
            return retval;
        } else {
            // split the sbend and insert a thinpole in between the halves
            ElmPtr sbptr1;
            ElmPtr sbptr2;
            elm->Split(0.5, sbptr1, sbptr2);
            std::vector < std::complex<double > > c_moments;
            for (int k = 0; k <= highest_order; ++k) {
                c_moments.push_back(std::complex<double >(bk[k], ak[k]));
            }

            retval.push_back(sbptr1);
            retval.push_back(
                    ElmPtr(
                            new ThinPole(
                                    (lattice_element.get_name() + "_poles").c_str(),
                                    brho * angle, c_moments)));
            retval.push_back(sbptr2);

            return retval;
        }
    } else {
        // combined function element
        CF_sbend* elm = new CF_sbend(lattice_element.get_name().c_str(),
                length, brho * angle / length, angle, e1, e2);
        // Does the CHEF default number of kicks match the Synergia default number of kicks?
        if (elm->numberOfKicks() != kicks) {
            elm->setNumberOfKicks(kicks);
        }
        if (tilt != 0.0) elm->setAlignment(aligner);
        double multipoleStrength = k1 * brho * length;
        if (multipoleStrength != 0.0) {
            elm->setQuadrupole(multipoleStrength);
        }
        multipoleStrength = k2 * brho * length / 2.0;
        if (multipoleStrength != 0.0) {
            elm->setSextupole(multipoleStrength);
        }
        multipoleStrength = k3 * brho * length / 6.0;
        if (multipoleStrength != 0.0) {
            elm->setOctupole(multipoleStrength);
        }
        ElmPtr elmP(elm);
        retval.push_back(elmP);
        return retval;
    }

}

template<class Archive>
    void
    Sbend_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Sbend_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Sbend_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Sbend_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Sbend_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Sbend_mad8_adaptor::~Sbend_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Sbend_mad8_adaptor)

Rbend_mad8_adaptor::Rbend_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("angle", 0.0);
    get_default_element().set_double_attribute("k1", 0.0);
    get_default_element().set_double_attribute("e1", 0.0);
    get_default_element().set_double_attribute("e2", 0.0);
    get_default_element().set_double_attribute("k2", 0.0);
    get_default_element().set_double_attribute("h1", 0.0);
    get_default_element().set_double_attribute("h2", 0.0);
    get_default_element().set_double_attribute("hgap", 0.0);
    get_default_element().set_double_attribute("fint", 0.0);
    get_default_element().set_double_attribute("k3", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
    // possible higher order multipole components
    get_default_element().set_double_attribute("kl", 0.0); // base strength/B-rho
    get_default_element().set_double_attribute("a1", 0.0); // skew quad
    get_default_element().set_double_attribute("a2", 0.0); // skew sextupole
    get_default_element().set_double_attribute("a3", 0.0); // skew octupole
    get_default_element().set_double_attribute("a4", 0.0); // skew decapole
    get_default_element().set_double_attribute("a5", 0.0); // skew dodecapole
    get_default_element().set_double_attribute("a6", 0.0); // skew tetradecapole
    get_default_element().set_double_attribute("a7", 0.0); // skew hexdecapole
    get_default_element().set_double_attribute("b1", 0.0); // quad
    get_default_element().set_double_attribute("b2", 0.0); // sextupole
    get_default_element().set_double_attribute("b3", 0.0); // octopole
    get_default_element().set_double_attribute("b4", 0.0); // decapole
    get_default_element().set_double_attribute("b5", 0.0); // dodecapole
    get_default_element().set_double_attribute("b6", 0.0); // tetradecapole
    get_default_element().set_double_attribute("b7", 0.0); // hexdecapole

}

void
Rbend_mad8_adaptor::set_defaults(Lattice_element & lattice_element)
{
    lattice_element.set_length_attribute_name("arclength");
    lattice_element.set_needs_internal_derive(true);
    Element_adaptor::set_defaults(lattice_element);
}

void
Rbend_mad8_adaptor::set_derived_attributes_internal(
        Lattice_element & lattice_element)
{
    double bend_angle = lattice_element.get_bend_angle();
    double bend_length = lattice_element.get_double_attribute("l");
    double arc_length = bend_angle * bend_length
            / (2 * std::sin(bend_angle / 2));
    lattice_element.set_double_attribute("arclength", arc_length);
}

Chef_elements
Rbend_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_double_attribute("l");
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");
    double k1 = lattice_element.get_double_attribute("k1");
    double k2 = lattice_element.get_double_attribute("k2");
    double k3 = lattice_element.get_double_attribute("k3");
    double tilt = lattice_element.get_double_attribute("tilt");
    bool simple = ((k1 == 0.0) && (k2 == 0.0) && (k3 == 0.0) && (tilt == 0.0));

    double ak[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently a0-a7
    double bk[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently b0-b7

    // thinpole strengths
    string a_attr_list[] = { "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7" };
    string b_attr_list[] = { "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7" };

    bool has_multipoles = false;
    int highest_order = 0;
    // find any possible multipole moments
    for (int moment = 0; moment < 8; ++moment) {
        if (lattice_element.has_double_attribute(a_attr_list[moment])) {
            ak[moment] = lattice_element.get_double_attribute(
                    a_attr_list[moment]);
            // no point in setting using multipole machinery ff the attributes are 0
            if (ak[moment] != 0.0) {
                has_multipoles = true;
                highest_order = moment;
            }
        }
        if (lattice_element.has_double_attribute(b_attr_list[moment])) {
            bk[moment] = lattice_element.get_double_attribute(
                    b_attr_list[moment]);
            if (bk[moment] != 0.0) {
                has_multipoles = true;
                highest_order = moment;
            }
        }
    }

    // mad8 implements rbends as parallel faced sbends.  Unless the
    // the string attribute true_rbend is "true", do the same.
    bool true_rbend = false;
    double arc_length = length;
    if (lattice_element.has_string_attribute("true_rbend") &&
        lattice_element.get_string_attribute("true_rbend") == "true") {
        true_rbend = true;
    } else {
        arc_length = length * angle/ (2 * std::sin(angle / 2));
    }

    // being not simple precludes using multipoles
    if (!simple && has_multipoles) {
        throw runtime_error("shouldn't use k1,k2 and multipoles in one rbend");
    }

    if (simple) {
        bmlnElmnt * bmelmnt;

        if ((0.0 == e1) && (0.0 == e2)) {
            if (true_rbend) {
                bmelmnt = new rbend(lattice_element.get_name().c_str(), length,
                                    brho * (2.0 * sin(0.5 * angle)) / length, angle);
                bmelmnt->setTag("RBEND");
            } else {
                bmelmnt = new sbend(lattice_element.get_name().c_str(), arc_length,
                                    brho * angle / arc_length, angle,
                                    angle/2.0, angle/2.0);
            }
        } else {
            bmelmnt = new rbend(lattice_element.get_name().c_str(), length,
                    brho * (2.0 * sin(0.5 * angle)) / length, angle, e1, e2);
            bmelmnt->setTag("RBEND");
        }
        ElmPtr elm(bmelmnt);
        // if there are no multipoles, I'm done.
        if (!has_multipoles) {
            retval.push_back(elm);
            return retval;
        } else {
            // split the sbend and insert a thinpole in between the halves
            ElmPtr rbptr1;
            ElmPtr rbptr2;
            elm->Split(0.5, rbptr1, rbptr2);

            std::vector < std::complex<double > > c_moments;
            for (int k = 0; k <= highest_order; ++k) {
                c_moments.push_back(std::complex<double >(bk[k], ak[k]));
            }

            retval.push_back(rbptr1);
            retval.push_back(
                    ElmPtr(
                            new ThinPole(
                                    (lattice_element.get_name() + "_poles").c_str(),
                                    brho * (2.0 * sin(0.5 * angle)),
                                    c_moments)));
            retval.push_back(rbptr2);

            return retval;
        }
    } else {
        // Not so simple
        alignmentData aligner;
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = tilt;

        bmlnElmnt* elm;
        if ((0.0 == e1) && (0.0 == e2)) {
            if (true_rbend) {
                elm = new CF_rbend(
                        lattice_element.get_name().c_str(), length,
                        brho * (2.0 * sin(0.5 * angle)) / length, angle);
            } else {
                elm = new CF_sbend(
                        lattice_element.get_name().c_str(), arc_length,
                        brho * angle / arc_length, angle,
                        angle/2.0, angle/2.0);
            }
        }
        else elm = new CF_rbend(lattice_element.get_name().c_str(), length,
                brho * (2.0 * sin(0.5 * angle)) / length, angle, e1, e2);

        elm->setTag("RBEND");
        if (tilt != 0.0) elm->setAlignment(aligner);

        double multipoleStrength = k1 * brho * length;
        if (multipoleStrength != 0.0) {
            if (true_rbend) {
                dynamic_cast<CF_rbend* >(elm)->setQuadrupole(multipoleStrength);
            } else {
                dynamic_cast<CF_sbend* >(elm)->setQuadrupole(multipoleStrength);
            }
        }
        multipoleStrength = k2 * brho * length / 2.0;
        if (multipoleStrength != 0.0) {
            if (true_rbend) {
                dynamic_cast<CF_rbend* >(elm)->setSextupole(multipoleStrength);
            } else {
                dynamic_cast<CF_sbend* >(elm)->setSextupole(multipoleStrength);
            }
        }
        multipoleStrength = k3 * brho * length / 6.0;
        if (multipoleStrength != 0.0) {
            if (true_rbend) {
                dynamic_cast<CF_rbend* >(elm)->setOctupole(multipoleStrength);
            } else {
                dynamic_cast<CF_sbend* >(elm)->setOctupole(multipoleStrength);
            }
        }

        ElmPtr elmP(elm);
        retval.push_back(elmP);
        return retval;
    }
}

template<class Archive>
    void
    Rbend_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Rbend_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rbend_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rbend_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rbend_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rbend_mad8_adaptor::~Rbend_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rbend_mad8_adaptor)

const char Quadrupole_mad8_adaptor::yoshida_propagator[] = "yoshida";
const char Quadrupole_mad8_adaptor::basic_propagator[] = "basic";

Quadrupole_mad8_adaptor::Quadrupole_mad8_adaptor()
{
    get_default_element().set_string_attribute("propagator_type", basic_propagator);
    get_default_element().set_double_attribute("yoshida_order", 2); // method is accurate to O[(kL)^(2*order+2)]
    get_default_element().set_double_attribute("yoshida_steps", 4);
    get_default_element().set_double_attribute("basic_kicks", 40);
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("k1", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
    get_default_element().set_double_attribute("hoffset", 0.0);
    get_default_element().set_double_attribute("voffset", 0.0);
    // possible higher order multipole components
    get_default_element().set_double_attribute("kl", 0.0); // base strength/B-rho
    get_default_element().set_double_attribute("a1", 0.0); // skew quad
    get_default_element().set_double_attribute("a2", 0.0); // skew sextupole
    get_default_element().set_double_attribute("a3", 0.0); // skew octupole
    get_default_element().set_double_attribute("a4", 0.0); // skew decapole
    get_default_element().set_double_attribute("a5", 0.0); // skew dodecapole
    get_default_element().set_double_attribute("a6", 0.0); // skew tetradecapole
    get_default_element().set_double_attribute("a7", 0.0); // skew hexdecapole
    get_default_element().set_double_attribute("b1", 0.0); // quad
    get_default_element().set_double_attribute("b2", 0.0); // sextupole
    get_default_element().set_double_attribute("b3", 0.0); // octopole
    get_default_element().set_double_attribute("b4", 0.0); // decapole
    get_default_element().set_double_attribute("b5", 0.0); // dodecapole
    get_default_element().set_double_attribute("b6", 0.0); // tetradecapole
    get_default_element().set_double_attribute("b7", 0.0); // hexdecapole
}

Chef_elements
Quadrupole_mad8_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{
    double qtilt;

    Chef_elements retval;

    alignmentData aligner;
    double length = lattice_element.get_double_attribute("l");

    double ak[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently a0-a7
    double bk[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently b0-b7

    // thinpole strengths
    string a_attr_list[] = { "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7" };
    string b_attr_list[] = { "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7" };

    double xoffset = lattice_element.get_double_attribute("hoffset");
    double yoffset = lattice_element.get_double_attribute("voffset");

    // a string attribute implies default pi/4
    if (lattice_element.has_string_attribute("tilt")) {
        qtilt = mconstants::pi / 4.0;
    } else {
        qtilt = lattice_element.get_double_attribute("tilt");
    }

    bool has_multipoles = false;
    int highest_order = 0;
    // find any possible multipole moments
    for (int moment = 0; moment < 8; ++moment) {
        if (lattice_element.has_double_attribute(a_attr_list[moment])) {
            ak[moment] = lattice_element.get_double_attribute(
                    a_attr_list[moment]);
            // no point in setting using multipole machinery ff the attributes are 0
            if (ak[moment] != 0.0) {
                has_multipoles = true;
                highest_order = moment;
            }
        }
        if (lattice_element.has_double_attribute(b_attr_list[moment])) {
            bk[moment] = lattice_element.get_double_attribute(
                    b_attr_list[moment]);
            if (bk[moment] != 0.0) {
                has_multipoles = true;
                highest_order = moment;
            }
        }
    }

    ElmPtr elm;
    if (length == 0.0) {
        elm.reset(new thinQuad(lattice_element.get_name().c_str(),
                brho * lattice_element.get_double_attribute("k1")));
    } else {
        if(lattice_element.get_string_attribute("propagator_type") == yoshida_propagator) {
            int steps = floor(lattice_element.get_double_attribute("yoshida_steps"));
            int order = floor(lattice_element.get_double_attribute("yoshida_order"));
            elm.reset(new quadrupole(lattice_element.get_name().c_str(), length,
                                        brho * lattice_element.get_double_attribute("k1")));
            quadrupole::PropagatorPtr yoshida_propagator(new YoshidaPropagator(order, steps));
            boost::dynamic_pointer_cast<quadrupole>(elm)->usePropagator(yoshida_propagator);
        } else if (lattice_element.get_string_attribute("propagator_type") == basic_propagator) {
            elm.reset(new quadrupole(lattice_element.get_name().c_str(), length,
                                        brho * lattice_element.get_double_attribute("k1")));
            boost::dynamic_pointer_cast<quadrupole>(elm)->setNumberOfKicks(floor(lattice_element.get_double_attribute("basic_kicks")));
        } else {
            throw std::runtime_error(
                        "Quadrupole_mad8_adaptor::get_chef_elements: bad propagator_type \"" +
                        lattice_element.get_string_attribute("propagator_type") + "\"");
        }
    }

    // using tilt and multipoles is a no-no
    if (has_multipoles && (qtilt != 0.0)) {
        throw runtime_error(
                "shouldn't use tilt and multipoles in same element");
    }

    bool needs_aligner;
    if ((qtilt != 0.0) || (xoffset != 0.0) || (yoffset != 0.0)) {
        needs_aligner = true;
        aligner.xOffset = xoffset;
        aligner.yOffset = yoffset;
        aligner.tilt = qtilt;
    } else {
        needs_aligner = false;
    }

    if (!has_multipoles) {
        if (needs_aligner) {
            elm->setAlignment(aligner);
        }
        retval.push_back(elm);
    } else {
        // split the quadrupole, insert thin multipole element in between halves
        std::vector < std::complex<double > > c_moments;
        for (int k = 0; k <= highest_order; ++k) {
            c_moments.push_back(std::complex<double >(bk[k], ak[k]));
        }

        double brkl = brho * length
                * lattice_element.get_double_attribute("k1");
        ElmPtr qptr1;
        ElmPtr qptr2;
        elm->Split(0.5, qptr1, qptr2);

        if (needs_aligner) {
            qptr1->setAlignment(aligner);
        }
        retval.push_back(qptr1);
        ElmPtr thinpoleptr(new ThinPole(
                (lattice_element.get_name() + "_poles").c_str(), brkl,
                c_moments));
        if (needs_aligner) {
            thinpoleptr->setAlignment(aligner);
        }
        retval.push_back(thinpoleptr);
        if (needs_aligner) {
            qptr2->setAlignment(aligner);
        }
        retval.push_back(qptr2);
    }

    return retval;
}

template<class Archive>
    void
    Quadrupole_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Quadrupole_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Quadrupole_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Quadrupole_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Quadrupole_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Quadrupole_mad8_adaptor::~Quadrupole_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Quadrupole_mad8_adaptor)

Sextupole_mad8_adaptor::Sextupole_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("k2", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Sextupole_mad8_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{
    Chef_elements retval;

    double sexlen = lattice_element.get_double_attribute("l");
    double sexk2 = lattice_element.get_double_attribute("k2");
    double sextilt = 0.0;

    alignmentData aligner;

    bmlnElmnt* bmln_elmnt;

    if (sexlen == 0.0) {
        bmln_elmnt = new thinSextupole(lattice_element.get_name().c_str(),
                brho * sexk2 / 2.0);
    } else {
        bmln_elmnt = new sextupole(lattice_element.get_name().c_str(), sexlen,
                brho * sexk2 / 2.0);
    }

    if (lattice_element.has_double_attribute("tilt")) {
        sextilt = lattice_element.get_double_attribute("tilt");
    } else if (lattice_element.has_string_attribute("tilt")) {
        // if this is a string, assume just tilt specified with no
        // value so use pi/6.
        sextilt = mconstants::pi / 6.0;
    }

    if (sextilt != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = sextilt;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Sextupole_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Sextupole_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Sextupole_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Sextupole_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Sextupole_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Sextupole_mad8_adaptor::~Sextupole_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Sextupole_mad8_adaptor)

Octupole_mad8_adaptor::Octupole_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("k3", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Octupole_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double octulen = lattice_element.get_double_attribute("l");
    double octuk2 = lattice_element.get_double_attribute("k3");
    double octutilt = 0.0;

    alignmentData aligner;

    bmlnElmnt* bmln_elmnt;

    if (octulen == 0.0) {
        bmln_elmnt = new thinOctupole(lattice_element.get_name().c_str(),
                brho * octuk2 / 6.0);
    } else {
        bmln_elmnt = new octupole(lattice_element.get_name().c_str(), octulen,
                brho * octuk2 / 6.0);
    }

    if (lattice_element.has_double_attribute("tilt")) {
        octutilt = lattice_element.get_double_attribute("tilt");
    } else if (lattice_element.has_string_attribute("tilt")) {
        // if this is a string, assume just tilt specified with no
        // value so use pi/8.
        octutilt = mconstants::pi / 8.0;
    }

    if (octutilt != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = octutilt;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Octupole_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Octupole_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Octupole_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Octupole_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Octupole_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Octupole_mad8_adaptor::~Octupole_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Octupole_mad8_adaptor)

Multipole_mad8_adaptor::Multipole_mad8_adaptor()
{
    get_default_element().set_double_attribute("k0l", 0.0);
    get_default_element().set_double_attribute("t0", 0.0);
    get_default_element().set_double_attribute("k1l", 0.0);
    get_default_element().set_double_attribute("t1", 0.0);
    get_default_element().set_double_attribute("k2l", 0.0);
    get_default_element().set_double_attribute("t2", 0.0);
    get_default_element().set_double_attribute("k3l", 0.0);
    get_default_element().set_double_attribute("t3", 0.0);
    get_default_element().set_double_attribute("k4l", 0.0);
    get_default_element().set_double_attribute("t4", 0.0);
    get_default_element().set_double_attribute("k5l", 0.0);
    get_default_element().set_double_attribute("t5", 0.0);
    get_default_element().set_double_attribute("k6l", 0.0);
    get_default_element().set_double_attribute("t6", 0.0);
    get_default_element().set_double_attribute("k7l", 0.0);
    get_default_element().set_double_attribute("t7", 0.0);
    get_default_element().set_double_attribute("k8l", 0.0);
    get_default_element().set_double_attribute("t8", 0.0);
    get_default_element().set_double_attribute("k9l", 0.0);
    get_default_element().set_double_attribute("t9", 0.0);
}

Chef_elements
Multipole_mad8_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{
    Chef_elements retval;

    // multipole strengths
    static string k_attr_list[] = { "k0l", "k1l", "k2l", "k3l", "k4l", "k5l",
            "k6l", "k7l", "k8l", "k9l" };
    // multipole tilts
    static string t_attr_list[] = { "t0", "t1", "t2", "t3", "t4", "t5", "t6",
            "t7", "t8", "t9" };

    double knl[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double tn[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    static double nfactorial[] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0,
            5040.0, 40320.0, 362880.0 };

    // loop through possible attributes
    for (int moment = 0; moment < 10; ++moment) {
        if (lattice_element.has_double_attribute(k_attr_list[moment])) {
            knl[moment] = lattice_element.get_double_attribute(
                    k_attr_list[moment]);
            // look for a tilt
            if (lattice_element.has_double_attribute(t_attr_list[moment])) {
                tn[moment] = lattice_element.get_double_attribute(
                        t_attr_list[moment]);
            } else if (lattice_element.has_string_attribute(
                    t_attr_list[moment])) {
                // any string value is equivalent to just giving the Tn keyword
                // get default tilt for that multipole order
                tn[moment] = mconstants::pi / (2 * moment + 2);
            }
        }
    }

    // assemble chef thinpoles for each specified knl multipole moment
    int multipole_count = 0;
    for (int moment = 0; moment < 10; ++moment) {

        if (knl[moment] != 0.0) {
            bmlnElmnt* bmln_elmnt = 0;
            alignmentData aligner;
            ++multipole_count;
            std::stringstream element_name(stringstream::out);

            element_name << lattice_element.get_name() << "_" << 2 * moment + 2
                    << "pole";
            bmln_elmnt = new ThinPole(element_name.str().c_str(),
                    brho * knl[moment] / nfactorial[moment], 2 * moment + 2);

            ElmPtr elm(bmln_elmnt);
            // set tilt if necessary
            if (tn[moment] != 0.0) {
                aligner.xOffset = 0.0;
                aligner.yOffset = 0.0;
                aligner.tilt = tn[moment];
                elm->setAlignment(aligner);
            }
            retval.push_back(elm);
        } else {
            std::stringstream element_name(stringstream::out);
            element_name << lattice_element.get_name() << "_" << 2 * moment + 2
                    << "pole_marker";
            ElmPtr elm = ElmPtr(new marker(element_name.str().c_str()));
            retval.push_back(elm);
        }
    }
    // csp: temporally or permanently disabled this part to avoid confusion.
    // put in a marker for this element
    //ElmPtr elm = ElmPtr(new marker(lattice_element.get_name().c_str()));
    //retval.push_back(elm);

    return retval;
}

template<class Archive>
    void
    Multipole_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Multipole_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Multipole_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Multipole_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Multipole_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Multipole_mad8_adaptor::~Multipole_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Multipole_mad8_adaptor)

//------------------------------------------
// thinpoles are an addon present only in CHEF
// they have length 0, normal and skew multipole moments
// they are specified by the kl factor of their enclosing dipole
//   or quadrupole, and b_k, and a_k coefficients relative to the
//   base element strength

Thinpole_mad8_adaptor::Thinpole_mad8_adaptor()
{
    get_default_element().set_double_attribute("kl", 0.0); // base strength/B-rho
    get_default_element().set_double_attribute("a1", 0.0); // skew quad
    get_default_element().set_double_attribute("a2", 0.0); // skew sextupole
    get_default_element().set_double_attribute("a3", 0.0); // skew octupole
    get_default_element().set_double_attribute("a4", 0.0); // skew decapole
    get_default_element().set_double_attribute("a5", 0.0); // skew dodecapole
    get_default_element().set_double_attribute("a6", 0.0); // skew tetradecapole
    get_default_element().set_double_attribute("a7", 0.0); // skew hexdecapole
    get_default_element().set_double_attribute("b1", 0.0); // quad
    get_default_element().set_double_attribute("b2", 0.0); // sextupole
    get_default_element().set_double_attribute("b3", 0.0); // octopole
    get_default_element().set_double_attribute("b4", 0.0); // decapole
    get_default_element().set_double_attribute("b5", 0.0); // dodecapole
    get_default_element().set_double_attribute("b6", 0.0); // tetradecapole
    get_default_element().set_double_attribute("b7", 0.0); // hexdecapole
}

Chef_elements
Thinpole_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    double ak[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently a0-a7
    double bk[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // currently b0-b7

    // thinpole strengths
    string a_attr_list[] = { "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7" };
    string b_attr_list[] = { "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7" };

    // loop through possible attributes
    for (int moment = 0; moment < 8; ++moment) {
        if (lattice_element.has_double_attribute(a_attr_list[moment])) {
            ak[moment] = lattice_element.get_double_attribute(
                    a_attr_list[moment]);
        }
        if (lattice_element.has_double_attribute(b_attr_list[moment])) {
            bk[moment] = lattice_element.get_double_attribute(
                    b_attr_list[moment]);
        }
    }

    double kl = lattice_element.get_double_attribute("kl");

    // assemble chef elements
    std::vector < std::complex<double > > c_moments;
    for (int k = 0; k < 8; ++k) {
        c_moments.push_back(std::complex<double >(bk[k], ak[k]));
    }

    retval.push_back(
            ElmPtr(
                    new ThinPole(lattice_element.get_name().c_str(), brho * kl,
                            c_moments)));

    // put in a marker for this element
    ElmPtr elm = ElmPtr(new marker(lattice_element.get_name().c_str()));
    retval.push_back(elm);

    return retval;
}

template<class Archive>
    void
    Thinpole_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Thinpole_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Thinpole_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Thinpole_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Thinpole_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Thinpole_mad8_adaptor::~Thinpole_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Thinpole_mad8_adaptor)

Solenoid_mad8_adaptor::Solenoid_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("ks", 0.0);
}

template<class Archive>
    void
    Solenoid_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Solenoid_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Solenoid_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Solenoid_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Solenoid_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Solenoid_mad8_adaptor::~Solenoid_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Solenoid_mad8_adaptor)

Hkicker_mad8_adaptor::Hkicker_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("kick", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Hkicker_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    alignmentData aligner;

    double kicklen = lattice_element.get_double_attribute("l");
    double kick = lattice_element.get_double_attribute("kick");
    double tilt = lattice_element.get_double_attribute("tilt");

    bmlnElmnt* bmln_elmnt = new hkick(lattice_element.get_name().c_str());
    if (kicklen > 0.0) {
        bmln_elmnt->setLength(kicklen);
    }
    if (kick != 0.0) {
        bmln_elmnt->setStrength(kick * brho);
    }
    if (tilt != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = tilt;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Hkicker_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Hkicker_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Hkicker_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Hkicker_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Hkicker_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Hkicker_mad8_adaptor::~Hkicker_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Hkicker_mad8_adaptor)

Vkicker_mad8_adaptor::Vkicker_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("kick", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Vkicker_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    alignmentData aligner;

    double kicklen = lattice_element.get_double_attribute("l");
    double kick = lattice_element.get_double_attribute("kick");
    double tilt = lattice_element.get_double_attribute("tilt");

    bmlnElmnt* bmln_elmnt = new vkick(lattice_element.get_name().c_str());
    if (kicklen > 0.0) {
        bmln_elmnt->setLength(kicklen);
    }
    if (kick != 0.0) {
        bmln_elmnt->setStrength(kick * brho);
    }
    if (tilt != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = tilt;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Vkicker_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Vkicker_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Vkicker_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Vkicker_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Vkicker_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Vkicker_mad8_adaptor::~Vkicker_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Vkicker_mad8_adaptor)

Kicker_mad8_adaptor::Kicker_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("hkick", 0.0);
    get_default_element().set_double_attribute("vkick", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Kicker_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    alignmentData aligner;

    double kicklen = lattice_element.get_double_attribute("l");
    double hkick = lattice_element.get_double_attribute("hkick");
    double vkick = lattice_element.get_double_attribute("vkick");
    double tilt = lattice_element.get_double_attribute("tilt");

    kick* bmln_elmnt = new kick(lattice_element.get_name().c_str());
    if (kicklen > 0.0) {
        bmln_elmnt->setLength(kicklen);
    }
    if (hkick != 0.0) {
        bmln_elmnt->setHorStrength(hkick * brho);
    }
    if (vkick != 0.0) {
        bmln_elmnt->setVerStrength(vkick * brho);
    }
    if (tilt > 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = tilt;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Kicker_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Kicker_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Kicker_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Kicker_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Kicker_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Kicker_mad8_adaptor::~Kicker_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Kicker_mad8_adaptor)

Rfcavity_mad8_adaptor::Rfcavity_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("volt", 0.0);
    get_default_element().set_double_attribute("lag", 0.0);
    get_default_element().set_double_attribute("harmon", 0.0);
    get_default_element().set_double_attribute("betrf", 0.0);
    get_default_element().set_double_attribute("pg", 0.0);
    get_default_element().set_double_attribute("shunt", 0.0);
    get_default_element().set_double_attribute("tfill", 0.0);
}

void
Rfcavity_mad8_adaptor::set_defaults(Lattice_element & lattice_element)
{
#if 0 // let CHEF set RF frequency
    lattice_element.set_needs_external_derive(true);
#endif // let CHEF set RF frequency
    Element_adaptor::set_defaults(lattice_element);
}

#if 0 // let CHEF set RF frequency
void
Rfcavity_mad8_adaptor::set_derived_attributes_external(Lattice_element &lattice_element,
		double lattice_length, double beta)
{
    if (lattice_element.has_double_attribute("harmon")
            && lattice_element.get_double_attribute("harmon") != 0.0) {
    	double h = lattice_element.get_double_attribute("harmon");
    	// freq in MHz to match what the input would be in a MAD file.
    	// I wish we didn't have to divide and multiply by 1.0e6.
    	double freq = 1.0e-6 * h * beta * pconstants::c/lattice_length;
    	lattice_element.set_double_attribute("freq", freq);
    }
}
#endif // let CHEF set RF frequency

Chef_elements
Rfcavity_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double freq = 0.0;

    if (lattice_element.has_double_attribute("freq")) {
    	freq = lattice_element.get_double_attribute("freq");
    }

    int harmonic_number = lattice_element.get_double_attribute("harmon");
    double q = 0;
	// Although mad8 does not support freq.  madx has freq in MHz.  We'll go
	// with the madx convention, converting to Hz for the CHEF constructor.
    // The CHEF rfcavity constructor takes voltage argument in eV, but
    // converts in to GeV for internal use. mad units for volt are MV.
    if (length == 0.0) {
        bmlnElmnt *bmln_elmnt;
        bmln_elmnt = new thinrfcavity(lattice_element.get_name().c_str(), freq*1.0e6,
                lattice_element.get_double_attribute("volt") * 1.0e6,
                lattice_element.get_double_attribute("lag")
                        * (2.0 * mconstants::pi), q,
                lattice_element.get_double_attribute("shunt"));
        static_cast<thinrfcavity *>(bmln_elmnt)->setHarmonicNumber(harmonic_number);
        ElmPtr elm(bmln_elmnt);
        retval.push_back(elm);
    } else {
        bmlnElmnt *pre_drift, *kick, *post_drift;
        pre_drift = new drift(
                (lattice_element.get_name() + "_predrift").c_str(),
                0.5 * length);
        kick = new thinrfcavity((lattice_element.get_name() + "_kick").c_str(),
                freq*1.0e6, lattice_element.get_double_attribute("volt") * 1.0e6,
                lattice_element.get_double_attribute("lag")
                        * (2.0 * mconstants::pi), q,
                lattice_element.get_double_attribute("shunt"));
        static_cast<thinrfcavity *>(kick)->setHarmonicNumber(harmonic_number);
        post_drift = new drift(
                (lattice_element.get_name() + "_postdrift").c_str(),
                0.5 * length);
        ElmPtr elm1(pre_drift), elm2(kick), elm3(post_drift);
        retval.push_back(elm1);
        retval.push_back(elm2);
        retval.push_back(elm3);
    }

    return retval;
}

template<class Archive>
    void
    Rfcavity_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Rfcavity_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rfcavity_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rfcavity_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rfcavity_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rfcavity_mad8_adaptor::~Rfcavity_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rfcavity_mad8_adaptor)

Elseparator_mad8_adaptor::Elseparator_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("e", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

template<class Archive>
    void
    Elseparator_mad8_adaptor::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Elseparator_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Elseparator_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Elseparator_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Elseparator_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Elseparator_mad8_adaptor::~Elseparator_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Elseparator_mad8_adaptor)

Hmonitor_mad8_adaptor::Hmonitor_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Hmonitor_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    bmlnElmnt* bmln_elmnt;

    double lenmon = lattice_element.get_double_attribute("l");
    bmln_elmnt = new hmonitor(lattice_element.get_name().c_str());
    if (lenmon > 0.0) {
        bmln_elmnt->setLength(lenmon);
    }

    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);

    return retval;
}

template<class Archive>
    void
    Hmonitor_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Hmonitor_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Hmonitor_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Hmonitor_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Hmonitor_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Hmonitor_mad8_adaptor::~Hmonitor_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Hmonitor_mad8_adaptor)

Vmonitor_mad8_adaptor::Vmonitor_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Vmonitor_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    bmlnElmnt* bmln_elmnt;

    double lenmon = lattice_element.get_double_attribute("l");
    bmln_elmnt = new vmonitor(lattice_element.get_name().c_str());
    if (lenmon > 0.0) {
        bmln_elmnt->setLength(lenmon);
    }

    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);

    return retval;
}

template<class Archive>
    void
    Vmonitor_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Vmonitor_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Vmonitor_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Vmonitor_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Vmonitor_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Vmonitor_mad8_adaptor::~Vmonitor_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Vmonitor_mad8_adaptor)

Monitor_mad8_adaptor::Monitor_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Monitor_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    bmlnElmnt* bmln_elmnt;

    double lenmon = lattice_element.get_double_attribute("l");
    bmln_elmnt = new monitor(lattice_element.get_name().c_str());
    if (lenmon > 0.0) {
        bmln_elmnt->setLength(lenmon);
    }

    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);

    return retval;
}

template<class Archive>
    void
    Monitor_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Monitor_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Monitor_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Monitor_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Monitor_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Monitor_mad8_adaptor::~Monitor_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Monitor_mad8_adaptor)

Instrument_mad8_adaptor::Instrument_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

// Use a drift to reserve space for instruments
Chef_elements
Instrument_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    bmlnElmnt* bmln_elmnt;

    double lenmon = lattice_element.get_double_attribute("l");
    bmln_elmnt = new drift(lattice_element.get_name().c_str(), lenmon);

    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);

    return retval;
}

template<class Archive>
    void
    Instrument_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Instrument_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Instrument_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Instrument_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Instrument_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Instrument_mad8_adaptor::~Instrument_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Instrument_mad8_adaptor)

Ecollimator_mad8_adaptor::Ecollimator_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("xsize", 0.0);
    get_default_element().set_double_attribute("ysize", 0.0);
}

Chef_elements
Ecollimator_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(
            new drift(lattice_element.get_name().c_str(),
                    lattice_element.get_length()));
    // add the Synergia aperture attributes?
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Ecollimator_mad8_adaptor::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Ecollimator_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Ecollimator_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Ecollimator_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Ecollimator_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Ecollimator_mad8_adaptor::~Ecollimator_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Ecollimator_mad8_adaptor)

Rcollimator_mad8_adaptor::Rcollimator_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("xsize", 0.0);
    get_default_element().set_double_attribute("ysize", 0.0);
}

Chef_elements
Rcollimator_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(
            new drift(lattice_element.get_name().c_str(),
                    lattice_element.get_length()));
    // add synergia aperture attributes?
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Rcollimator_mad8_adaptor::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Rcollimator_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rcollimator_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rcollimator_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rcollimator_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rcollimator_mad8_adaptor::~Rcollimator_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rcollimator_mad8_adaptor)

Septum_mad8_adaptor::Septum_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("voltage", 0.0);
    get_default_element().set_double_attribute("gap", 0.0);
    get_default_element().set_double_attribute("wire_position", 0.0);
    get_default_element().set_double_attribute("wire_width", 0.0);
    get_default_element().set_double_attribute("positive_strength", 0.0);
    get_default_element().set_double_attribute("negative_strength", 0.0);
}

Chef_elements
Septum_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_double_attribute("l");
    double voltage = lattice_element.get_double_attribute("voltage");
    double gap = lattice_element.get_double_attribute("gap");
    double wire_position = lattice_element.get_double_attribute(
            "wire_position");
    double wire_width = lattice_element.get_double_attribute("wire_width");
    double positive_strength = lattice_element.get_double_attribute(
            "positive_strength");
    double negative_strength = lattice_element.get_double_attribute(
            "negative_strength");

    if (length == 0.0) {
        ThinSeptumPtr es_septum(
                new thinSeptum(lattice_element.get_name().c_str(),
                        positive_strength, negative_strength, wire_position));
        es_septum->setStrengths(positive_strength, negative_strength);
        es_septum->setWire(wire_position);
        es_septum->setWireWidth(wire_width);
        es_septum->setGap(gap);
        ElmPtr elm(es_septum);
        //TODO: below doesn't work. why?
        //bmlnElmnt * bmln_elmnt;
        //bmln_elmnt = new thinSeptum(lattice_element.get_name().c_str(),
        //        positive_strength, negative_strength, wire_position);
        //bmln_elmnt->setStrengths(positive_strength, negative_strength);
        //bmln_elmnt->setWire(wire_position);
        //bmln_elmnt->setWireWidth(wire_width);
        //bmln_elmnt->setGap(gap);
        //ElmPtr elm(bmln_elmnt);
        retval.push_back(elm);
    } else {
        SeptumPtr es_septum(
                new septum(lattice_element.get_name().c_str(), length, voltage,
                        gap, wire_position, wire_width));
        ElmPtr elm(es_septum);
        retval.push_back(elm);
    }
    return retval;
}

template<class Archive>
    void
    Septum_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Septum_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Septum_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Septum_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Septum_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Septum_mad8_adaptor::~Septum_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Septum_mad8_adaptor)

Lambertson_mad8_adaptor::Lambertson_mad8_adaptor()
{
}

Chef_elements
Lambertson_mad8_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{
    Chef_elements retval;

    bmlnElmnt * bmln_elmnt;
    bmln_elmnt = new thinLamb(lattice_element.get_name().c_str());
    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Lambertson_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Lambertson_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lambertson_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lambertson_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lambertson_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Lambertson_mad8_adaptor::~Lambertson_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Lambertson_mad8_adaptor)


Srot_mad8_adaptor::Srot_mad8_adaptor()
{
}


Chef_elements
Srot_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    double angle = lattice_element.get_double_attribute("angle");

    ElmPtr elm(
            new srot(lattice_element.get_name().c_str(), angle));

    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Srot_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Srot_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Srot_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Srot_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Srot_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);



Srot_mad8_adaptor::~Srot_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Srot_mad8_adaptor)

Elens_mad8_adaptor::Elens_mad8_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("current", 0.0);
    get_default_element().set_double_attribute("eenergy", 0.0);
    get_default_element().set_double_attribute("radius", 0.0);
    get_default_element().set_double_attribute("gaussian", 0.0);
    get_default_element().set_double_attribute("uniform", 0.0);
}

Chef_elements
Elens_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double lenslen = lattice_element.get_double_attribute("l");
    double current = lattice_element.get_double_attribute("current");
    double eenergy = lattice_element.get_double_attribute("eenergy")*0.001;  // convert eenergy from MV to GV
    double radius = lattice_element.get_double_attribute("radius");
    bool gaussian = !(lattice_element.get_double_attribute("gaussian") == 0.0);
    bool uniform = !(lattice_element.get_double_attribute("uniform") == 0.0);


    if (!(uniform || gaussian)) {
        throw std::runtime_error("elens must set either gaussian or uniform attribute");
    }
    if (gaussian && uniform) {
        throw std::runtime_error("elens must not set both gaussian and uniform attributes");
    }

    elens::e_profile_t prof(elens::undefined);
    if (gaussian) {
        prof = elens::gaussian;
    } else if (uniform) {
        prof = elens::uniform;
    }
    ElensPtr elensptr(new elens(lattice_element.get_name().c_str(), lenslen, current, eenergy, radius, prof));

    retval.push_back(elensptr);
    return retval;
}

template<class Archive>
    void
    Elens_mad8_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Elens_mad8_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Elens_mad8_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Elens_mad8_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Elens_mad8_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Elens_mad8_adaptor::~Elens_mad8_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Elens_mad8_adaptor)
