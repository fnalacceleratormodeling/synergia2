#include "element_adaptor.h"

#include <beamline/beamline_elements.h>
#include "synergia/foundation/math_constants.h"

Element_adaptor::Element_adaptor()
{
}

void
Element_adaptor::set_double_default(Lattice_element & lattice_element,
        std::string const& name, double value)
{
    if (!lattice_element.has_double_attribute(name)) {
        lattice_element.set_double_attribute(name, value);
    }
}

void
Element_adaptor::set_string_default(Lattice_element & lattice_element,
        std::string const& name, std::string const& value)
{
    if (!lattice_element.has_string_attribute(name)) {
        lattice_element.set_string_attribute(name, value);
    }
}

void
Element_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
}

Chef_elements
Element_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    throw(runtime_error("Element_adaptor: " + lattice_element.get_type()
            + " not handled"));
}

Element_adaptor::~Element_adaptor()
{
}

Element_adaptor_map::Element_adaptor_map()
{
    boost::shared_ptr<Marker_mad8_adaptor > marker_mad8_adaptor(
            new Marker_mad8_adaptor);
    adaptor_map["marker"] = marker_mad8_adaptor;

    boost::shared_ptr<Drift_mad8_adaptor > drift_mad8_adaptor(
            new Drift_mad8_adaptor);
    adaptor_map["drift"] = drift_mad8_adaptor;

    boost::shared_ptr<Sbend_mad8_adaptor > sbend_mad8_adaptor(
            new Sbend_mad8_adaptor);
    adaptor_map["sbend"] = sbend_mad8_adaptor;

    boost::shared_ptr<Rbend_mad8_adaptor > rbend_mad8_adaptor(
            new Rbend_mad8_adaptor);
    adaptor_map["rbend"] = rbend_mad8_adaptor;

    boost::shared_ptr<Quadrupole_mad8_adaptor > quadrupole_mad8_adaptor(
            new Quadrupole_mad8_adaptor);
    adaptor_map["quadrupole"] = quadrupole_mad8_adaptor;

    boost::shared_ptr<Sextupole_mad8_adaptor > sextupole_mad8_adaptor(
            new Sextupole_mad8_adaptor);
    adaptor_map["sextupole"] = sextupole_mad8_adaptor;

    boost::shared_ptr<Octupole_mad8_adaptor > octupole_mad8_adaptor(
            new Octupole_mad8_adaptor);
    adaptor_map["octupole"] = octupole_mad8_adaptor;

    boost::shared_ptr<Multipole_mad8_adaptor > multipole_mad8_adaptor(
            new Multipole_mad8_adaptor);
    adaptor_map["multipole"] = multipole_mad8_adaptor;

    boost::shared_ptr<Solenoid_mad8_adaptor > solenoid_mad8_adaptor(
            new Solenoid_mad8_adaptor);
    adaptor_map["solenoid"] = solenoid_mad8_adaptor;

    boost::shared_ptr<Hkicker_mad8_adaptor > hkicker_mad8_adaptor(
            new Hkicker_mad8_adaptor);
    adaptor_map["hkicker"] = hkicker_mad8_adaptor;

    boost::shared_ptr<Vkicker_mad8_adaptor > vkicker_mad8_adaptor(
            new Vkicker_mad8_adaptor);
    adaptor_map["vkicker"] = vkicker_mad8_adaptor;

    boost::shared_ptr<Kicker_mad8_adaptor > kicker_mad8_adaptor(
            new Kicker_mad8_adaptor);
    adaptor_map["kicker"] = kicker_mad8_adaptor;

    boost::shared_ptr<Rfcavity_mad8_adaptor > rfcavity_mad8_adaptor(
            new Rfcavity_mad8_adaptor);
    adaptor_map["rfcavity"] = rfcavity_mad8_adaptor;

    boost::shared_ptr<Elseparator_mad8_adaptor > elseparator_mad8_adaptor(
            new Elseparator_mad8_adaptor);
    adaptor_map["elseparator"] = elseparator_mad8_adaptor;

    boost::shared_ptr<Hmonitor_mad8_adaptor > hmonitor_mad8_adaptor(
            new Hmonitor_mad8_adaptor);
    adaptor_map["hmonitor"] = hmonitor_mad8_adaptor;

    boost::shared_ptr<Vmonitor_mad8_adaptor > vmonitor_mad8_adaptor(
            new Vmonitor_mad8_adaptor);
    adaptor_map["vmonitor"] = vmonitor_mad8_adaptor;

    boost::shared_ptr<Monitor_mad8_adaptor > monitor_mad8_adaptor(
            new Monitor_mad8_adaptor);
    adaptor_map["monitor"] = monitor_mad8_adaptor;

    boost::shared_ptr<Instrument_mad8_adaptor > instrument_mad8_adaptor(
            new Instrument_mad8_adaptor);
    adaptor_map["instrument"] = instrument_mad8_adaptor;

    boost::shared_ptr<Ecollimator_mad8_adaptor > ecollimator_mad8_adaptor(
            new Ecollimator_mad8_adaptor);
    adaptor_map["ecollimator"] = ecollimator_mad8_adaptor;

    boost::shared_ptr<Rcollimator_mad8_adaptor > rcollimator_mad8_adaptor(
            new Rcollimator_mad8_adaptor);
    adaptor_map["rcollimator"] = rcollimator_mad8_adaptor;
}

void
Element_adaptor_map::set_adaptor(std::string const& name,
        Element_adaptor_sptr element_adaptor_sptr)
{
    adaptor_map[name] = element_adaptor_sptr;
}

bool
Element_adaptor_map::has_adaptor(std::string const& name) const
{
    return (adaptor_map.count(name) > 0);
}

Element_adaptor_sptr
Element_adaptor_map::get_adaptor(std::string const& name) const
{
    std::map<std::string, Element_adaptor_sptr >::const_iterator iter =
            adaptor_map.find(name);
    return iter->second;
}

std::list<std::string >
Element_adaptor_map::get_adaptor_names() const
{
    std::list<std::string > retval;
    for (std::map<std::string, Element_adaptor_sptr >::const_iterator it =
            adaptor_map.begin(); it != adaptor_map.end(); ++it) {
        retval.push_back(it->first);
    }
    return retval;
}

Element_adaptor_map::~Element_adaptor_map()
{
}

Marker_mad8_adaptor::Marker_mad8_adaptor()
{
}

void
Marker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
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

Marker_mad8_adaptor::~Marker_mad8_adaptor()
{
}

Drift_mad8_adaptor::Drift_mad8_adaptor()
{
}

void
Drift_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
}

Chef_elements
Drift_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;
    ElmPtr elm(new drift(lattice_element.get_name().c_str(),
            lattice_element.get_length()));
    retval.push_back(elm);
    return retval;
}

Drift_mad8_adaptor::~Drift_mad8_adaptor()
{
}

Sbend_mad8_adaptor::Sbend_mad8_adaptor()
{
}

void
Sbend_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "angle", 0.0);
    set_double_default(lattice_element, "k1", 0.0);
    set_double_default(lattice_element, "e1", 0.0);
    set_double_default(lattice_element, "e2", 0.0);
    set_double_default(lattice_element, "k2", 0.0);
    set_double_default(lattice_element, "h1", 0.0);
    set_double_default(lattice_element, "h2", 0.0);
    set_double_default(lattice_element, "hgap", 0.0);
    set_double_default(lattice_element, "fint", 0.0);
    set_double_default(lattice_element, "k3", 0.0);
    if (!lattice_element.has_double_attribute("tilt")
            && !lattice_element.has_string_attribute("tilt")) {
        lattice_element.set_double_attribute("tilt", 0.0);
    }
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

    if ((lattice_element.get_double_attribute("k1") != 0.0)
            || (lattice_element.get_double_attribute("k2") != 0.0)
            || (lattice_element.get_double_attribute("k3") != 0.0)
            || (lattice_element.get_double_attribute("tilt") != 0.0)
            || (lattice_element.get_double_attribute("h1") != 0.0)
            || (lattice_element.get_double_attribute("h2") != 0.0)
            || (lattice_element.get_double_attribute("hgap") != 0.0)
            || (lattice_element.get_double_attribute("fint") != 0.0)) {
        throw(runtime_error(
                "lattice_element_to_chef_sbend: non-zero element(s) of something not handled"));
    }

    ElmPtr elm(new sbend(lattice_element.get_name().c_str(), length, brho
            * angle / length, angle, e1, e2));
    retval.push_back(elm);
    return retval;
}

Sbend_mad8_adaptor::~Sbend_mad8_adaptor()
{
}

Rbend_mad8_adaptor::Rbend_mad8_adaptor()
{
}

void
Rbend_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "angle", 0.0);
    set_double_default(lattice_element, "k1", 0.0);
    set_double_default(lattice_element, "e1", 0.0);
    set_double_default(lattice_element, "e2", 0.0);
    set_double_default(lattice_element, "k2", 0.0);
    set_double_default(lattice_element, "h1", 0.0);
    set_double_default(lattice_element, "h2", 0.0);
    set_double_default(lattice_element, "hgap", 0.0);
    set_double_default(lattice_element, "fint", 0.0);
    set_double_default(lattice_element, "k3", 0.0);
    if (!lattice_element.has_double_attribute("tilt")
            && !lattice_element.has_string_attribute("tilt")) {
        lattice_element.set_double_attribute("tilt", 0.0);
    }
}

Chef_elements
Rbend_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double angle = lattice_element.get_double_attribute("angle");
    double e1 = lattice_element.get_double_attribute("e1");
    double e2 = lattice_element.get_double_attribute("e2");

    if ((lattice_element.get_double_attribute("k1") != 0.0)
            || (lattice_element.get_double_attribute("k2") != 0.0)
            || (lattice_element.get_double_attribute("k3") != 0.0)
            || (lattice_element.get_double_attribute("tilt") != 0.0)
            || (lattice_element.get_double_attribute("h1") != 0.0)
            || (lattice_element.get_double_attribute("h2") != 0.0)
            || (lattice_element.get_double_attribute("hgap") != 0.0)
            || (lattice_element.get_double_attribute("fint") != 0.0)) {
        throw(runtime_error(
                "lattice_element_to_chef_sbend: non-zero element(s) of something not handled"));
    }

    ElmPtr elm(new rbend(lattice_element.get_name().c_str(), length, brho
            * (2.0 * sin(0.5 * angle)) / length, angle, e1, e2));
    retval.push_back(elm);
    return retval;
}

Rbend_mad8_adaptor::~Rbend_mad8_adaptor()
{
}

Quadrupole_mad8_adaptor::Quadrupole_mad8_adaptor()
{
}

void
Quadrupole_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "k1", 0.0);
    if (!lattice_element.has_double_attribute("tilt")
            && !lattice_element.has_string_attribute("tilt")) {
        lattice_element.set_double_attribute("tilt", 0.0);
    }

}

Chef_elements
Quadrupole_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    double qtilt;
    Chef_elements retval;

    alignmentData aligner;
    double length = lattice_element.get_double_attribute("l");

    // a string attribute implies default pi/4
    if (lattice_element.has_string_attribute("tilt")) {
        qtilt = M_PI / 4.0;
    } else {
        qtilt = lattice_element.get_double_attribute("tilt");
    }

    bmlnElmnt* bmln_elmnt;
    if (length == 0.0) {
        bmln_elmnt = new thinQuad(lattice_element.get_name().c_str(), brho
                * lattice_element.get_double_attribute("k1"));
    } else {
        bmln_elmnt = new quadrupole(lattice_element.get_name().c_str(), length,
                brho * lattice_element.get_double_attribute("k1"));
    }

    if (qtilt != 0.0) {
        aligner.xOffset = 0.0;
        aligner.yOffset = 0.0;
        aligner.tilt = qtilt;
        bmln_elmnt->setAlignment(aligner);
    }

    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

Quadrupole_mad8_adaptor::~Quadrupole_mad8_adaptor()
{
}

Sextupole_mad8_adaptor::Sextupole_mad8_adaptor()
{
}

void
Sextupole_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "k2", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Chef_elements
Sextupole_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double sexlen = lattice_element.get_double_attribute("l");
    double sexk2 = lattice_element.get_double_attribute("k2");
    double sextilt;

    alignmentData aligner;

    bmlnElmnt* bmln_elmnt;

    if (sexlen == 0.0) {
        bmln_elmnt = new thinSextupole(lattice_element.get_name().c_str(), brho
                * sexk2 / 2.0);
    } else {
        bmln_elmnt = new sextupole(lattice_element.get_name().c_str(), sexlen,
                brho * sexk2 / 2.0);
    }

    if (lattice_element.has_double_attribute("tilt")) {
        sextilt = lattice_element.get_double_attribute("tilt");
    } else if (lattice_element.has_string_attribute("tilt")) {
        // if this is a string, assume just tilt specified with no
        // value so use pi/6.
        sextilt = M_PI / 6.0;
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

Sextupole_mad8_adaptor::~Sextupole_mad8_adaptor()
{
}

Octupole_mad8_adaptor::Octupole_mad8_adaptor()
{
}

void
Octupole_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "k3", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Chef_elements
Octupole_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double octulen = lattice_element.get_double_attribute("l");
    double octuk2 = lattice_element.get_double_attribute("k3");
    double octutilt;

    alignmentData aligner;

    bmlnElmnt* bmln_elmnt;

    if (octulen == 0.0) {
        bmln_elmnt = new thinOctupole(lattice_element.get_name().c_str(), brho
                * octuk2 / 6.0);
    } else {
        bmln_elmnt = new octupole(lattice_element.get_name().c_str(), octulen,
                brho * octuk2 / 6.0);
    }

    if (lattice_element.has_double_attribute("tilt")) {
        octutilt = lattice_element.get_double_attribute("tilt");
    } else if (lattice_element.has_string_attribute("tilt")) {
        // if this is a string, assume just tilt specified with no
        // value so use pi/8.
        octutilt = M_PI / 8.0;
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
Octupole_mad8_adaptor::~Octupole_mad8_adaptor()
{
}

Multipole_mad8_adaptor::Multipole_mad8_adaptor()
{
}

void
Multipole_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "k0l", 0.0);
    set_double_default(lattice_element, "t0", 0.0);
    set_double_default(lattice_element, "k1l", 0.0);
    set_double_default(lattice_element, "t1", 0.0);
    set_double_default(lattice_element, "k2l", 0.0);
    set_double_default(lattice_element, "t2", 0.0);
    set_double_default(lattice_element, "k3l", 0.0);
    set_double_default(lattice_element, "t3", 0.0);
    set_double_default(lattice_element, "k4l", 0.0);
    set_double_default(lattice_element, "t4", 0.0);
    set_double_default(lattice_element, "k5l", 0.0);
    set_double_default(lattice_element, "t5", 0.0);
    set_double_default(lattice_element, "k6l", 0.0);
    set_double_default(lattice_element, "t6", 0.0);
    set_double_default(lattice_element, "k7l", 0.0);
    set_double_default(lattice_element, "t7", 0.0);
    set_double_default(lattice_element, "k8l", 0.0);
    set_double_default(lattice_element, "t8", 0.0);
    set_double_default(lattice_element, "k9l", 0.0);
    set_double_default(lattice_element, "t9", 0.0);
}

Chef_elements
Multipole_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    alignmentData aligner;

    // multipole strengths
    string k_attr_list[] = { "k0l", "k1l", "k2l", "k3l", "k4l", "k5l", "k6l",
            "k7l", "k8l", "k9l" };
    // multipole tilts
    string t_attr_list[] = { "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7",
            "t8", "t9" };

    double knl[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double tn[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // loop through possible attributes
    for (int moment = 0; moment < 10; ++moment) {
        if (lattice_element.has_double_attribute(k_attr_list[moment])) {
            knl[moment] = lattice_element.get_double_attribute(
                    k_attr_list[moment]);
            // look for a tilt
            if (lattice_element.has_double_attribute(t_attr_list[moment])) {
                tn[moment] = lattice_element.get_double_attribute(
                        t_attr_list[moment]);
            } else if (lattice_element.has_string_attribute(t_attr_list[moment])) {
                // any string value is equivalent to just giving the Tn keyword
                // get default tilt for that multipole order
                tn[moment] = M_PI / (2 * moment + 2);
            }
        }
    }

    // assemble chef elements
    int multipole_count = 0;
    for (int moment = 0; moment < 10; ++moment) {
        if (knl[moment] != 0.0) {
            bmlnElmnt* bmln_elmnt;
            ++multipole_count;
            switch (moment) {

            case 0:
                bmln_elmnt = new thin2pole("", brho * knl[moment]);
                break;
            case 1:
                bmln_elmnt = new thinQuad("", brho * knl[moment]);
                break;
            case 2:
                bmln_elmnt = new thinSextupole("", brho * knl[moment] / 2.0);
                break;
            case 3:
                bmln_elmnt = new thinOctupole("", brho * knl[moment] / 6.0);
                break;
            case 4:
                bmln_elmnt = new thinDecapole("", brho * knl[moment] / 24.0);
                break;
            }

            ElmPtr elm(bmln_elmnt);
            // set tilt if necessary
            if (tn[moment] != 0.0) {
                aligner.xOffset = 0.0;
                aligner.yOffset = 0.0;
                aligner.tilt = tn[moment];
                elm->setAlignment(aligner);
            }
            retval.push_back(elm);
        }
    }

    // put in a marker for this element
    ElmPtr elm = ElmPtr(new marker(lattice_element.get_name().c_str()));
    ;
    retval.push_back(elm);
    return retval;
}

Multipole_mad8_adaptor::~Multipole_mad8_adaptor()
{
}

Solenoid_mad8_adaptor::Solenoid_mad8_adaptor()
{
}

void
Solenoid_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "ks", 0.0);
}

Solenoid_mad8_adaptor::~Solenoid_mad8_adaptor()
{
}

Hkicker_mad8_adaptor::Hkicker_mad8_adaptor()
{
}

void
Hkicker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "kick", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
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

Hkicker_mad8_adaptor::~Hkicker_mad8_adaptor()
{
}

Vkicker_mad8_adaptor::Vkicker_mad8_adaptor()
{
}

void
Vkicker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "kick", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
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

Vkicker_mad8_adaptor::~Vkicker_mad8_adaptor()
{
}

Kicker_mad8_adaptor::Kicker_mad8_adaptor()
{
}

void
Kicker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "hkick", 0.0);
    set_double_default(lattice_element, "vkick", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
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

Kicker_mad8_adaptor::~Kicker_mad8_adaptor()
{
}

Rfcavity_mad8_adaptor::Rfcavity_mad8_adaptor()
{
}

void
Rfcavity_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "volt", 0.0);
    set_double_default(lattice_element, "lag", 0.0);
    set_double_default(lattice_element, "harmon", 0.0);
    set_double_default(lattice_element, "betrf", 0.0);
    set_double_default(lattice_element, "pg", 0.0);
    set_double_default(lattice_element, "shunt", 0.0);
    set_double_default(lattice_element, "tfill", 0.0);
}

Chef_elements
Rfcavity_mad8_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    Chef_elements retval;

    double length = lattice_element.get_length();
    double freq = 0;
    if (lattice_element.has_double_attribute("freq")) {
        freq = lattice_element.get_double_attribute("freq");
    } else {
        if (lattice_element.has_double_attribute("harmon")
                && lattice_element.get_double_attribute("harmon") != 0.0) {
            std::cout
                    << "jfa: rfcavity could figure out frequency from harmonic number, but doesn't. FIXME!\n";
        }
    }
    double q = 0;
    bmlnElmnt* bmln_elmnt;
    if (length == 0.0) {
        bmln_elmnt = new thinrfcavity(lattice_element.get_name().c_str(), freq,
                lattice_element.get_double_attribute("volt") * 1.0e6,
                lattice_element.get_double_attribute("lag") * (2.0
                        * mconstants::pi), q,
                lattice_element.get_double_attribute("shunt"));
    } else {
        bmln_elmnt = new rfcavity(lattice_element.get_name().c_str(), length,
                freq, lattice_element.get_double_attribute("volt") * 1.0e6,
                lattice_element.get_double_attribute("lag") * (2.0
                        * mconstants::pi), q,
                lattice_element.get_double_attribute("shunt"));
    }

    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}

Rfcavity_mad8_adaptor::~Rfcavity_mad8_adaptor()
{
}

Elseparator_mad8_adaptor::Elseparator_mad8_adaptor()
{
}

void
Elseparator_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "e", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Elseparator_mad8_adaptor::~Elseparator_mad8_adaptor()
{
}

Hmonitor_mad8_adaptor::Hmonitor_mad8_adaptor()
{
}

void
Hmonitor_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
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

Hmonitor_mad8_adaptor::~Hmonitor_mad8_adaptor()
{
}

Vmonitor_mad8_adaptor::Vmonitor_mad8_adaptor()
{
}

void
Vmonitor_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
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

Vmonitor_mad8_adaptor::~Vmonitor_mad8_adaptor()
{
}

Monitor_mad8_adaptor::Monitor_mad8_adaptor()
{
}

void
Monitor_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
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

Monitor_mad8_adaptor::~Monitor_mad8_adaptor()
{
}

Instrument_mad8_adaptor::Instrument_mad8_adaptor()
{
}

void
Instrument_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
}

Instrument_mad8_adaptor::~Instrument_mad8_adaptor()
{
}

Ecollimator_mad8_adaptor::Ecollimator_mad8_adaptor()
{
}

void
Ecollimator_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "xsize", 0.0);
    set_double_default(lattice_element, "ysize", 0.0);
}

Ecollimator_mad8_adaptor::~Ecollimator_mad8_adaptor()
{
}

Rcollimator_mad8_adaptor::Rcollimator_mad8_adaptor()
{
}

void
Rcollimator_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "xsize", 0.0);
    set_double_default(lattice_element, "ysize", 0.0);
}

Rcollimator_mad8_adaptor::~Rcollimator_mad8_adaptor()
{
}

