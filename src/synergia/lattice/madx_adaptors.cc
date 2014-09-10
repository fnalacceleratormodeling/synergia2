#include <stdexcept>
#include <sstream>
#include "madx_adaptors.h"
#include "synergia/foundation/math_constants.h"

#include <complex>

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
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

Marker_madx_adaptor::Marker_madx_adaptor()
{
}

Chef_elements
Marker_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(new marker(lattice_element.get_name().c_str()));
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Marker_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Marker_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Marker_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Marker_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Marker_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Marker_madx_adaptor::~Marker_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Marker_madx_adaptor)

Drift_madx_adaptor::Drift_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Drift_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Drift_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Drift_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Drift_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Drift_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Drift_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Drift_madx_adaptor::~Drift_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Drift_madx_adaptor)

Sbend_madx_adaptor::Sbend_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("angle", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
    get_default_element().set_double_attribute("k1", 0.0);
    get_default_element().set_double_attribute("e1", 0.0);
    get_default_element().set_double_attribute("e2", 0.0);
    get_default_element().set_double_attribute("fint", 0.0);
    get_default_element().set_double_attribute("fintx", 0.0);
    get_default_element().set_double_attribute("hgap", 0.0);
    get_default_element().set_double_attribute("k2", 0.0);
    get_default_element().set_double_attribute("h1", 0.0);
    get_default_element().set_double_attribute("h2", 0.0);
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
Sbend_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");
    double k1 = lattice_element.get_double_attribute("k1");
    double k2 = lattice_element.get_double_attribute("k2");
    double tilt = lattice_element.get_double_attribute("tilt");

    bool simple = ((k1 == 0.0) && (k2 == 0.0));

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
    	// std::cout << "sbend element: " << lattice_element.get_name() << ", angle: " << angle << ", e1: " << e1 << ", e2: " << e2 << ", k1: " << k1 << ", k2: " << k2 << std::endl;
    	ElmPtr elm;
    	// if angle is 0.0, this reduces to a drift, but only if the edge angles are also 0
    	if (angle == 0.0 && e1 == 0.0 && e2 == 0.0) {
    		elm = ElmPtr(
    				new drift((lattice_element.get_name()+"_conv2drft").c_str(), length));
    	} else {
    		elm = ElmPtr(
    				new sbend(lattice_element.get_name().c_str(), length,
    						brho * angle / length, angle, e1, e2));
    	}
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
        bmlnElmnt* elm = new CF_sbend(lattice_element.get_name().c_str(),
                length, brho * angle / length, angle, e1, e2);
        if (tilt != 0.0) elm->setAlignment(aligner);
        double multipoleStrength = k1 * brho * length;
        if (multipoleStrength != 0.0) {
            dynamic_cast<CF_sbend* >(elm)->setQuadrupole(multipoleStrength);
        }
        multipoleStrength = k2 * brho * length / 2.0;
        if (multipoleStrength != 0.0) {
            dynamic_cast<CF_sbend* >(elm)->setSextupole(multipoleStrength);
        }
        ElmPtr elmP(elm);
        retval.push_back(elmP);
        return retval;
    }

}

template<class Archive>
    void
    Sbend_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Sbend_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Sbend_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Sbend_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Sbend_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Sbend_madx_adaptor::~Sbend_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Sbend_madx_adaptor)

Rbend_madx_adaptor::Rbend_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("angle", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
    // jfa: add_angle not yet handled
    get_default_element().set_double_attribute("k1", 0.0);
    get_default_element().set_double_attribute("e1", 0.0);
    get_default_element().set_double_attribute("e2", 0.0);
    get_default_element().set_double_attribute("fint", 0.0);
    get_default_element().set_double_attribute("fintx", 0.0);
    get_default_element().set_double_attribute("hgap", 0.0);
    get_default_element().set_double_attribute("k2", 0.0);
    get_default_element().set_double_attribute("h1", 0.0);
    get_default_element().set_double_attribute("h2", 0.0);
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
Rbend_madx_adaptor::set_defaults(Lattice_element & lattice_element)
{
    lattice_element.set_length_attribute_name("arclength");
    lattice_element.set_needs_internal_derive(true);
    Element_adaptor::set_defaults(lattice_element);
}

void
Rbend_madx_adaptor::set_derived_attributes_internal(
        Lattice_element & lattice_element)
{
    double bend_angle = lattice_element.get_bend_angle();
    double bend_length = lattice_element.get_double_attribute("l");
    double arc_length = bend_angle * bend_length
            / (2 * std::sin(bend_angle / 2));
    lattice_element.set_double_attribute("arclength", arc_length);
}

Chef_elements
Rbend_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_double_attribute("l");
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");
    double k1 = lattice_element.get_double_attribute("k1");
    double k2 = lattice_element.get_double_attribute("k2");
    double tilt = lattice_element.get_double_attribute("tilt");
    bool simple = ((k1 == 0.0) && (k2 == 0.0) && (tilt == 0.0));

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

    // being not simple precludes using multipoles
    if (!simple && has_multipoles) {
        throw runtime_error("shouldn't use k1,k2 and multipoles in one rbend");
    }

    if (simple) {
        bmlnElmnt * bmelmnt;

        if ((0.0 == e1) && (0.0 == e2)) {
            bmelmnt = new rbend(lattice_element.get_name().c_str(), length,
                    brho * (2.0 * sin(0.5 * angle)) / length, angle);
            bmelmnt->setTag("RBEND");
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
        // std::cout << "rbend not so simple: " << lattice_element.get_name() << std::endl; std::cout.flush();
        // std::cout << "length: " << length << ", strength: " << brho * (2.0 * sin(0.5 * angle)) / length << std::endl; std::cout.flush();
        // Not so simple
        alignmentData aligner;
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = tilt;

        bmlnElmnt* elm;
        if ((0.0 == e1) && (0.0 == e2)) elm = new CF_rbend(
                lattice_element.get_name().c_str(), length,
                brho * (2.0 * sin(0.5 * angle)) / length, angle);
        else elm = new CF_rbend(lattice_element.get_name().c_str(), length,
                brho * (2.0 * sin(0.5 * angle)) / length, angle, e1, e2);

        elm->setTag("RBEND");
        if (tilt != 0.0) elm->setAlignment(aligner);

        double multipoleStrength = k1 * brho * length;
        if (multipoleStrength != 0.0) {
            dynamic_cast<CF_rbend* >(elm)->setQuadrupole(multipoleStrength);
        }
        multipoleStrength = k2 * brho * length / 2.0;
        if (multipoleStrength != 0.0) {
            dynamic_cast<CF_rbend* >(elm)->setSextupole(multipoleStrength);
        }

        ElmPtr elmP(elm);
        retval.push_back(elmP);
        return retval;
    }
}

template<class Archive>
    void
    Rbend_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Rbend_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rbend_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rbend_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rbend_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rbend_madx_adaptor::~Rbend_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rbend_madx_adaptor)

const char Quadrupole_madx_adaptor::yoshida_propagator[] = "yoshida";
const char Quadrupole_madx_adaptor::basic_propagator[] = "basic";

Quadrupole_madx_adaptor::Quadrupole_madx_adaptor()
{
    get_default_element().set_string_attribute("propagator_type", yoshida_propagator);
    get_default_element().set_double_attribute("yoshida_order", 2); // method is accurate to O[(kL)^(2*order+2)]
    get_default_element().set_double_attribute("yoshida_steps", 4);
    get_default_element().set_double_attribute("basic_kicks", 40);
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("k1", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
    get_default_element().set_double_attribute("k1s", 0.0);
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
Quadrupole_madx_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{

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

    double qtilt = lattice_element.get_double_attribute("tilt");

    // ck1 is the complex strength of the quadrupole k1 - i*k1s (negative because of MAD-X definition)
    std::complex<double> ck1(lattice_element.get_double_attribute("k1"),
                             -lattice_element.get_double_attribute("k1s"));
    // rot_angle is phase of ck1/(order+1) (order = 1 for quadrupole)
    // for pure skew quad, phase is -pi/2, order=1 so rotation = -pi/4
    // so rot_angle is the rotation needed to turn a normal quad
    // into one with k1 and k1s components.
     double rot_angle = std::arg(ck1)/2.0 + qtilt;

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

    bmlnElmnt* bmln_elmnt;
    if (length == 0.0) {
        bmln_elmnt = new thinQuad(lattice_element.get_name().c_str(),
                                  brho * std::abs(ck1));
    } else {
        if(lattice_element.get_string_attribute("propagator_type") == yoshida_propagator) {
            int steps = floor(lattice_element.get_double_attribute("yoshida_steps"));
            int order = floor(lattice_element.get_double_attribute("yoshida_order"));
            bmln_elmnt = new quadrupole(lattice_element.get_name().c_str(), length,
                                        brho * std::abs(ck1));
            quadrupole::PropagatorPtr yoshida_propagator(new YoshidaPropagator(order, steps));
            dynamic_cast<quadrupole*>(bmln_elmnt)->usePropagator(yoshida_propagator);
        } else if (lattice_element.get_string_attribute("propagator_type") == basic_propagator) {
            bmln_elmnt = new quadrupole(lattice_element.get_name().c_str(), length,
                                        brho * std::abs(ck1));
            dynamic_cast<quadrupole*>(bmln_elmnt)->setNumberOfKicks(floor(lattice_element.get_double_attribute("basic_kicks")));
        } else {
            throw std::runtime_error(
                        "Quadrupole_madx_adaptor::get_chef_elements: bad propagator_type \"" +
                        lattice_element.get_string_attribute("propagator_type") + "\"");
        }
    }

    bool needs_aligner;
    if ((rot_angle != 0.0) || (xoffset != 0.0) || (yoffset != 0.0)) {
        needs_aligner = true;
        aligner.xOffset = xoffset;
        aligner.yOffset = yoffset;
        aligner.tilt = rot_angle;
    } else {
        needs_aligner = false;
    }

    if (!has_multipoles) {
        if (needs_aligner) {
            bmln_elmnt->setAlignment(aligner);
        }
        ElmPtr elm(bmln_elmnt);
        retval.push_back(elm);
     } else {
        // split the quadrupole, insert thin multipole element in between halves
        std::vector < std::complex<double > > c_moments;
        for (int k = 0; k <= highest_order; ++k) {
            c_moments.push_back(std::complex<double >(bk[k], ak[k]));
        }
        // if explicit qtilt attribute set, rotation accomplished by multiplying
        // b + i*a coefficients by exp(-i (n+1)*theta ).  The multipoles are rotated by
        // the tilt angle, but the magnet body is rotated by normal+skew+tilt angle
        for (int k=0; k <= highest_order; ++k) {
            c_moments[k] *= std::complex<double>(std::cos((k+1)*qtilt), -std::sin((k+1)*qtilt));
        }

        double brkl = brho * length
                * lattice_element.get_double_attribute("k1");
        ElmPtr qptr1;
        ElmPtr qptr2;
        bmln_elmnt->Split(0.5, qptr1, qptr2);

        if (needs_aligner) {
            qptr1->setAlignment(aligner);
        }
        retval.push_back(qptr1);
        bmlnElmnt* thinpoleptr = new ThinPole(
                (lattice_element.get_name() + "_poles").c_str(), brkl,
                c_moments);
        if (needs_aligner) {
            thinpoleptr->setAlignment(aligner);
        }
        retval.push_back(ElmPtr(thinpoleptr));
        if (needs_aligner) {
            qptr2->setAlignment(aligner);
        }
        retval.push_back(qptr2);
    }

    return retval;
}

template<class Archive>
    void
    Quadrupole_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Quadrupole_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Quadrupole_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Quadrupole_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Quadrupole_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Quadrupole_madx_adaptor::~Quadrupole_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Quadrupole_madx_adaptor)

Sextupole_madx_adaptor::Sextupole_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("k2", 0.0);
    get_default_element().set_double_attribute("k2s", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Sextupole_madx_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{

    Chef_elements retval;

    double sexlen = lattice_element.get_double_attribute("l");
    alignmentData aligner;

    bmlnElmnt* bmln_elmnt;

    double sextilt = lattice_element.get_double_attribute("tilt");

    // ck2 is the complex strength of the sextupole k2 - i*k2s (negative because of MAD-X definition)
    std::complex<double> ck2(lattice_element.get_double_attribute("k2"),
                             -lattice_element.get_double_attribute("k2s"));
    // rot_angle is phase of ck2/(order+1) (order=2 for sextupole)
    // for pure skew sext, phase is -pi/2, order=2 so rotation = -pi/6
    // so rot_angle is the rotation needed to turn a normal sextupole
    // into one with k2 and k2s components.
     double rot_angle = std::arg(ck2)/3.0 + sextilt;

    if (sexlen == 0.0) {
        bmln_elmnt = new thinSextupole(lattice_element.get_name().c_str(),
                brho * std::abs(ck2) / 2.0);
    } else {
        bmln_elmnt = new sextupole(lattice_element.get_name().c_str(), sexlen,
                brho * std::abs(ck2) / 2.0);
    }

    if (rot_angle != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = rot_angle;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Sextupole_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Sextupole_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Sextupole_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Sextupole_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Sextupole_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Sextupole_madx_adaptor::~Sextupole_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Sextupole_madx_adaptor)

Octupole_madx_adaptor::Octupole_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("k3", 0.0);
    get_default_element().set_double_attribute("k3s", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Octupole_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{

    Chef_elements retval;

    double octulen = lattice_element.get_double_attribute("l");

    double octutilt = lattice_element.get_double_attribute("tilt");

    alignmentData aligner;

    bmlnElmnt* bmln_elmnt;

    // ck3 is the complex strength of the octupole k3 - i*k3s (negative because of MAD-X definition)
    std::complex<double> ck3(lattice_element.get_double_attribute("k3"),
                             -lattice_element.get_double_attribute("k3s"));
    // rot_angle is phase of ck3/(order+1) (order=3 for octupole)
    // for pure skew oct, phase is -pi/2, order=3 so rotation = -pi/8
    // so rot_angle is the rotation needed to turn a normal octupole
    // into one with k3 and k3s components.
     double rot_angle = std::arg(ck3)/4.0 + octutilt;

    if (octulen == 0.0) {
        bmln_elmnt = new thinOctupole(lattice_element.get_name().c_str(),
                brho * std::abs(ck3) / 6.0);
    } else {
        bmln_elmnt = new octupole(lattice_element.get_name().c_str(), octulen,
                brho * std::abs(ck3) / 6.0);
    }

    if (rot_angle != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = rot_angle;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm = ElmPtr(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Octupole_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Octupole_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Octupole_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Octupole_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Octupole_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Octupole_madx_adaptor::~Octupole_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Octupole_madx_adaptor)

Multipole_madx_adaptor::Multipole_madx_adaptor()
{
    get_default_element().set_double_attribute("tilt", 0.0);
    get_default_element().set_vector_attribute("knl", std::vector<double >(0));
    get_default_element().set_vector_attribute("ksl", std::vector<double >(0));
}

Chef_elements
Multipole_madx_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
{
    Chef_elements retval;

    double tilt = lattice_element.get_double_attribute("tilt");

    std::vector<double > knl(lattice_element.get_vector_attribute("knl"));
    double nfactorial = 1.0;
    for (int moment = 0; moment < knl.size(); ++moment) {
        std::stringstream element_name(stringstream::out);
        element_name << lattice_element.get_name() << "_" << 2 * moment + 2
                << "pole";
        ElmPtr elm(
                new ThinPole(element_name.str().c_str(),
                        brho * knl[moment] / nfactorial, 2 * moment + 2));
        if (tilt != 0.0) {
            alignmentData aligner;
            aligner.xOffset = 0.0;
            aligner.yOffset = 0.0;
            aligner.tilt = tilt;
            elm->setAlignment(aligner);
        }
        retval.push_back(elm);
        nfactorial *= (moment + 1);
    }

    std::vector<double > ksl(lattice_element.get_vector_attribute("ksl"));
    nfactorial = 1.0;
    for (int moment = 0; moment < ksl.size(); ++moment) {
        std::stringstream element_name(stringstream::out);
        element_name << lattice_element.get_name() << "_skew" << 2 * moment + 2
                << "pole";
        ElmPtr elm(
                new ThinPole(element_name.str().c_str(),
                        brho * ksl[moment] / nfactorial, 2 * moment + 2));
        alignmentData aligner;
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = mconstants::pi / (2 * moment + 2) + tilt;
        elm->setAlignment(aligner);
        retval.push_back(elm);
        nfactorial *= (moment + 1);
    }
    return retval;
}

template<class Archive>
    void
    Multipole_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Multipole_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Multipole_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Multipole_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Multipole_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Multipole_madx_adaptor::~Multipole_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Multipole_madx_adaptor)

//------------------------------------------
// thinpoles are an addon present only in CHEF
// they have length 0, normal and skew multipole moments
// they are specified by the kl factor of their enclosing dipole
//   or quadrupole, and b_k, and a_k coefficients relative to the
//   base element strength

Thinpole_madx_adaptor::Thinpole_madx_adaptor()
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
Thinpole_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Thinpole_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Thinpole_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Thinpole_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Thinpole_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Thinpole_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Thinpole_madx_adaptor::~Thinpole_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Thinpole_madx_adaptor)

Solenoid_madx_adaptor::Solenoid_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("ks", 0.0);
    get_default_element().set_double_attribute("ksi", 0.0);
}

Chef_elements
Solenoid_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
                                        double brho)
{
    Chef_elements retval;
    double length = lattice_element.get_double_attribute("l");
    double ks = lattice_element.get_double_attribute("ks");

    if (length == 0.0) {
        throw 
        std::runtime_error("Solenoid_madx_adaptor: zero-length solenoids not yet handled");
    }
    bmlnElmnt * generic_elm;
    if (ks == 0.0) {
        generic_elm = new drift(lattice_element.get_name().c_str(), length);
    } else {
        generic_elm = new Solenoid(lattice_element.get_name().c_str(), length, ks);
    }
    ElmPtr elm(generic_elm);
    retval.push_back(elm);
    return retval;
}

template<class Archive>
    void
    Solenoid_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Solenoid_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Solenoid_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Solenoid_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Solenoid_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Solenoid_madx_adaptor::~Solenoid_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Solenoid_madx_adaptor)

Hkicker_madx_adaptor::Hkicker_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("kick", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Hkicker_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Hkicker_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Hkicker_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Hkicker_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Hkicker_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Hkicker_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Hkicker_madx_adaptor::~Hkicker_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Hkicker_madx_adaptor)

Vkicker_madx_adaptor::Vkicker_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("kick", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Vkicker_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Vkicker_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Vkicker_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Vkicker_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Vkicker_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Vkicker_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Vkicker_madx_adaptor::~Vkicker_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Vkicker_madx_adaptor)

Kicker_madx_adaptor::Kicker_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("hkick", 0.0);
    get_default_element().set_double_attribute("vkick", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

Chef_elements
Kicker_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Kicker_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Kicker_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Kicker_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Kicker_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Kicker_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Kicker_madx_adaptor::~Kicker_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Kicker_madx_adaptor)

Rfcavity_madx_adaptor::Rfcavity_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("volt", 0.0);
    get_default_element().set_double_attribute("lag", 0.0);
    get_default_element().set_double_attribute("harmon", 0.0);
    get_default_element().set_double_attribute("shunt", 0.0);
}

void
Rfcavity_madx_adaptor::set_defaults(Lattice_element & lattice_element)
{
#if 0 // let CHEF set RF frequency
    lattice_element.set_needs_external_derive(true);
#endif // let CHEF set RF frequency
    Element_adaptor::set_defaults(lattice_element);
}

#if 0 // let CHEF set the RF frequency
void
Rfcavity_madx_adaptor::set_derived_attributes_external(Lattice_element &lattice_element,
		double lattice_length, double beta)
{
    if (lattice_element.has_double_attribute("harmon")
            && lattice_element.get_double_attribute("harmon") != 0.0) {
    	double h = lattice_element.get_double_attribute("harmon");

    	double freq = 1.0e-6 * h * beta * pconstants::c/lattice_length;
    	lattice_element.set_double_attribute("freq", freq);
    }
}
#endif // let CHEF set the RF frequency

Chef_elements
Rfcavity_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double freq = 0;
    if (lattice_element.has_double_attribute("freq")) {
        freq = lattice_element.get_double_attribute("freq");
    }
    int harmonic_number = lattice_element.get_double_attribute("harmon");
    double q = 0;
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
    Rfcavity_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Rfcavity_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rfcavity_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rfcavity_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rfcavity_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rfcavity_madx_adaptor::~Rfcavity_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rfcavity_madx_adaptor)

Elseparator_madx_adaptor::Elseparator_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("e", 0.0);
    get_default_element().set_double_attribute("tilt", 0.0);
}

template<class Archive>
    void
    Elseparator_madx_adaptor::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Elseparator_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Elseparator_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Elseparator_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Elseparator_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Elseparator_madx_adaptor::~Elseparator_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Elseparator_madx_adaptor)

Hmonitor_madx_adaptor::Hmonitor_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Hmonitor_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Hmonitor_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Hmonitor_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Hmonitor_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Hmonitor_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Hmonitor_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Hmonitor_madx_adaptor::~Hmonitor_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Hmonitor_madx_adaptor)

Vmonitor_madx_adaptor::Vmonitor_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Vmonitor_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Vmonitor_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Vmonitor_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Vmonitor_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Vmonitor_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Vmonitor_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Vmonitor_madx_adaptor::~Vmonitor_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Vmonitor_madx_adaptor)

Monitor_madx_adaptor::Monitor_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Monitor_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Monitor_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Monitor_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Monitor_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Monitor_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Monitor_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Monitor_madx_adaptor::~Monitor_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Monitor_madx_adaptor)

Instrument_madx_adaptor::Instrument_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
}

Chef_elements
Instrument_madx_adaptor::get_chef_elements(
        Lattice_element const& lattice_element, double brho)
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
    Instrument_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Instrument_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Instrument_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Instrument_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Instrument_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Instrument_madx_adaptor::~Instrument_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Instrument_madx_adaptor)

Ecollimator_madx_adaptor::Ecollimator_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("xsize", 0.0);
    get_default_element().set_double_attribute("ysize", 0.0);
}

template<class Archive>
    void
    Ecollimator_madx_adaptor::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Ecollimator_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Ecollimator_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Ecollimator_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Ecollimator_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Ecollimator_madx_adaptor::~Ecollimator_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Ecollimator_madx_adaptor)

Rcollimator_madx_adaptor::Rcollimator_madx_adaptor()
{
    get_default_element().set_double_attribute("l", 0.0);
    get_default_element().set_double_attribute("xsize", 0.0);
    get_default_element().set_double_attribute("ysize", 0.0);
}

template<class Archive>
    void
    Rcollimator_madx_adaptor::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Rcollimator_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rcollimator_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rcollimator_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rcollimator_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rcollimator_madx_adaptor::~Rcollimator_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rcollimator_madx_adaptor)

Septum_madx_adaptor::Septum_madx_adaptor()
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
Septum_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Septum_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Septum_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Septum_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Septum_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Septum_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Septum_madx_adaptor::~Septum_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Septum_madx_adaptor)

Lambertson_madx_adaptor::Lambertson_madx_adaptor()
{
}

Chef_elements
Lambertson_madx_adaptor::get_chef_elements(
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
    Lambertson_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Lambertson_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lambertson_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lambertson_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lambertson_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Lambertson_madx_adaptor::~Lambertson_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Lambertson_madx_adaptor)


Srot_madx_adaptor::Srot_madx_adaptor()
{
}


Chef_elements
Srot_madx_adaptor::get_chef_elements(Lattice_element const& lattice_element,
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
    Srot_madx_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor);
    }

template
void
Srot_madx_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Srot_madx_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Srot_madx_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Srot_madx_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);



Srot_madx_adaptor::~Srot_madx_adaptor()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Srot_madx_adaptor)
