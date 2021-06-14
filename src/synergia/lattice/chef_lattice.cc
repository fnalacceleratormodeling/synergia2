#include "chef_lattice.h"
#include "chef_utils.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/math_constants.h"
#include "chef_lattice_section.h"

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
#include <basic_toolkit/PhysicsConstants.h>
#include <physics_toolkit/DriftConverter.h>
#include <beamline/RefRegVisitor.h>
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

#include <stdexcept>
#include <string>
#include <cstring>

template<class Archive>
    void
    Chef_lattice::Begin_end::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(begin);
        ar & BOOST_SERIALIZATION_NVP(end);
    }

template
void
Chef_lattice::Begin_end::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Chef_lattice::Begin_end::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Chef_lattice::Begin_end::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Chef_lattice::Begin_end::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

const char Chef_lattice::internal_marker_name[] =
        "_synergia_lattice_element_marker";

void
Chef_lattice::construct_beamline()
{
    BmlPtr unpolished_beamline_sptr(new beamline());
    for (Lattice_elements::const_iterator latt_it =
            lattice_sptr->get_elements().begin(); latt_it
            != lattice_sptr->get_elements().end(); ++latt_it) {
        std::string type((*latt_it)->get_type());
        if (!lattice_sptr->get_element_adaptor_map().has_adaptor(type)) {
            throw(runtime_error("Chef_lattice: " + type + " not handled"));
        } else {
            Chef_elements
                    celms =
                            lattice_sptr->get_element_adaptor_map().get_adaptor(type)->get_chef_elements(
                                    *(*latt_it), brho);
            for (Chef_elements::const_iterator cel_it = celms.begin(); cel_it
                    != celms.end(); ++cel_it) {
                unpolished_beamline_sptr->append(*cel_it);
            }
            unpolished_beamline_sptr->append(lattice_element_marker);
        }
    }
    beamline_sptr = polish_beamline(unpolished_beamline_sptr);
    extract_element_map();
}

struct strengthData {
    thinrfcavity *address;
    double strength;
};

void
Chef_lattice::register_beamline(BmlPtr beamline_sptr)
{

    Particle testpart(reference_particle_to_chef_particle(lattice_sptr->get_reference_particle()));
    // std::cout << "Registering beamline with particle state: " << testpart.State() << std::endl;

    RefRegVisitor registrar(testpart);

    beamline_sptr->accept(registrar);

}

// is_a_bend is a free function
bool
is_a_bend(ElmPtr ce)
{
    return (strcmp(ce->Type(), "sbend") == 0 ||
            strcmp(ce->Type(), "rbend") == 0 ||
            strcmp(ce->Type(), "CF_sbend") == 0 ||
            strcmp(ce->Type(), "CF_rbend") == 0);
}

// find next bend starting from idx.  Return the index of the bend or -1.
int
find_next_bend(std::vector<beamline::iterator> const& biters, int idx, bool wraparound)
{
    int startidx = idx;
    bool do_wrap = wraparound && (startidx == 0);
    while (idx != biters.size()) {
        if (is_a_bend(*biters[idx])) {
            return idx;
        } else {
            ++idx;
        }
    }
    if (!do_wrap) {
        return -1; // didn't find it and we're not wrapping around
    }
    idx = 0;
    while (idx != startidx) {
        if (is_a_bend(*biters[idx])) {
            return idx;
        } else {
            ++idx;
        }
    }
    return -1; // didn't find it after wrap-around
}

// find next Slot that goes with the current bend starting from idx
int find_next_Slot(std::vector<beamline::iterator> const & biters, int idx, std::string findname, bool wraparound) {
    int startidx = idx;
    bool do_wrap = wraparound && (startidx == 0);
    while (idx != biters.size()) {
        if (strcmp((*biters[idx])->Type(), "Slot") == 0 &&
            (*biters[idx])->Name() == findname) {
            return idx;
        } else {
            ++idx;
        }
    }
    if (!do_wrap) {
        return -1; // we didn't find it and we're not wrapping around
    }
    idx = 0;
    while (idx != startidx) {
        if (strcmp((*biters[idx])->Type(), "Slot") == 0 &&
            (*biters[idx])->Name() == findname) {
            return idx;
        } else {
            ++idx;
        }
    }
    return -1; //we didn't find it even after wrapping around
}

// find previous Slot that goes with the current bend starting from index idx
int
find_prev_Slot(std::vector<beamline::iterator> const & biters, int idx, std::string findname, bool wraparound) {
    int startidx = idx;
    bool do_wrap = wraparound && (startidx == 0);
    while (idx > 0) {
        if (strcmp((*biters[idx])->Type(), "Slot") == 0 &&
            (*biters[idx])->Name() == findname) {
            return idx;
        } else {
            --idx;
        }
    }
    if (!do_wrap) {
        return -1; // we didn't find it and we're not wrapping around
    }
    idx = biters.size()-1;
    while (idx != startidx) {
        if (strcmp((*biters[idx])->Type(), "Slot") == 0 &&
            (*biters[idx])->Name() == findname) {
            return idx;
        } else {
            --idx;
        }
    }
    return -1; // we didn't find it even after wrapping around
}

void print_neighborhood(std::vector <beamline::iterator> const& biters, int idx, int n=3);

void print_neighborhood(std::vector<beamline::iterator> const & biters, int idx, int n)
{
    int startidx = (idx-n >= 0) ? (idx-n) : 0;
    int endidx = (idx+n < biters.size()) ? idx+n : biters.size();
    for (int p=startidx; p<endidx; ++p) {
        if (p==idx) {
            std::cout << "---> ";
        } else {
            std::cout << "     ";
        }
        std::cout << p << ": " << (*biters[p])->Name() << "(";
        std::cout << (*biters[p])->Type() << ")" << std::endl;
    }
    return;
}

BmlPtr
Chef_lattice::polish_beamline(BmlPtr beamline_sptr)
{
    DriftConverter drift_converter;
    BmlPtr converted_beamline_sptr;
    if (beamline_sptr->countHowManyDeeply() < 3) {
        converted_beamline_sptr = beamline_sptr;
    } else {
        converted_beamline_sptr = drift_converter.convert(*beamline_sptr);
    }

    // reorder beamline slots so chef elements
    //    slot | synergia_marker | magnet | slot | synergia_marker
    // becomes
    //    synergia_marker | slot | magnet | slot | synergia_marker

    // make a vector of all the beamline iterators for convenient indexing
    std::vector<beamline::iterator> biters;
    for (auto it=converted_beamline_sptr->begin();
         it!=converted_beamline_sptr->end(); ++it) {
        biters.push_back(it);
    }

    // std::cout << "polish_beamline:on entry" << std::endl;
    // print_neighborhood(biters, 0, 20);
    // std::cout << std::endl;

    // search for bends
    int cur_elem = 0;
    while (cur_elem < biters.size()) {
        int next_bend = find_next_bend(biters, cur_elem, false);
        if (next_bend < 0) {
            break; // no next bend, I'm done
        } else {
            cur_elem = next_bend;
        }
        
        // std::cout << "found bend at idx: " << cur_elem << std::endl;
        // print_neighborhood(biters, cur_elem, 4);
        // check for marker just before?
        if ((cur_elem > 0) && ((*biters[cur_elem-1])->Name() == lattice_element_marker->Name())) {
            int prev_slot = find_prev_Slot(biters, cur_elem, (*biters[cur_elem])->Name()+"_inSlot", false);
            // Move the prev_slot to just after the marker if it exists
            if (prev_slot >=0 ) {
                // std::cout << "Found previous slot at idx: " << prev_slot << std::endl;
                beamline::iterator inslotit = biters[prev_slot];
                for (int p=prev_slot+1; p != cur_elem; ++p) {
                    biters[p-1] = biters[p];
                }
                biters[cur_elem-1] = inslotit;
                // std::cout << "after rearrange: " << std::endl;
                // print_neighborhood(biters, cur_elem, 4);
            }
        }

        // check for marker afterwards
        if ((cur_elem < biters.size()-1) && ((*biters[cur_elem+1])->Name() == lattice_element_marker->Name())) {
            int next_slot = find_next_Slot(biters, cur_elem, (*biters[cur_elem])->Name()+"_outSlot", false);
            // Move the next_slot to just before the marker if it exists
            if (next_slot >= 0) {
                // std::cout << "Found next slot at idx: " << next_slot << std::endl;
                beamline::iterator outslotit = biters[next_slot];
                for (int p=next_slot-1; p != cur_elem; --p) {
                    biters[p+1] = biters[p];
                }
                biters[cur_elem+1] = outslotit;
                // std::cout << "after rearrange: " << std::endl;
                // print_neighborhood(biters, cur_elem, 4);
            }
        }
        // std::cout << std::endl;

        ++cur_elem;
    }

    BmlPtr slots_reordered_sptr(new beamline("slots-reordered-beamline"));
    for (auto& bit:biters) {
        slots_reordered_sptr->append(*bit);
    }
    register_beamline(slots_reordered_sptr);
    return slots_reordered_sptr;
}

void
Chef_lattice::extract_element_map()
{
    Lattice_elements::const_iterator le_it = lattice_sptr->get_elements().begin();
    Begin_end begin_end;
    begin_end.begin = 0;
    int size = beamline_sptr->countHowMany();
    beamline_iterators.resize(size + 1);
    int index = 0;
    for (beamline::iterator it = beamline_sptr->begin(); it
            != beamline_sptr->end(); ++it) {
        beamline_iterators.at(index) = it;
        if ((*it)->Name() == lattice_element_marker->Name()) {
            begin_end.end = index;
            element_map[le_it->get()] = begin_end;
            begin_end.begin = index + 1;
            ++le_it;
        }
        ++index;
    }
    beamline_iterators.at(index) = beamline_sptr->end();
}

void
Chef_lattice::extract_element_slice_map(Lattice_element_slices const& slices)
{
    Lattice_element_slices::const_iterator slice_it = slices.begin();
    Begin_end begin_end;
    begin_end.begin = 0;
    int size = sliced_beamline_sptr->countHowMany();
    sliced_beamline_iterators.resize(size + 1);
    sliced_beamline_const_iterators.resize(size + 1);
    int index = 0;
    beamline::const_iterator cit = sliced_beamline_sptr->begin();
    for (beamline::iterator it = sliced_beamline_sptr->begin(); it
            != sliced_beamline_sptr->end(); ++it) {
        sliced_beamline_iterators.at(index) = it;
        sliced_beamline_const_iterators.at(index) = cit;
        if ((*it)->Name() == lattice_element_marker->Name()) {
            begin_end.end = index;
            element_slice_map[slice_it->get()] = begin_end;
            begin_end.begin = index + 1;
            ++slice_it;
        }
        ++index;
        ++cit;
    }
    sliced_beamline_iterators.at(index) = sliced_beamline_sptr->end();
    sliced_beamline_const_iterators.at(index) = sliced_beamline_sptr->end();
}

void
Chef_lattice::construct()
{
    lattice_sptr->complete_attributes();
    sliced_beamline_sptr = BmlPtr(new beamline("sliced"));
    have_sliced_beamline_ = false;
    if (!lattice_sptr->has_reference_particle()) {
        throw(std::runtime_error(
                "Chef_lattice: requires a reference particle in Lattice"));
    }
    brho = lattice_sptr->get_reference_particle().get_momentum()
            / PH_CNV_brho_to_p;

    construct_beamline();
}

Chef_lattice::Chef_lattice(Lattice_sptr lattice_sptr) :
        lattice_sptr(lattice_sptr), beamline_sptr(), lattice_element_marker(
                new marker(internal_marker_name))
{
    construct();
}

double
Chef_lattice::get_brho() const
{
    return brho;
}

Chef_lattice::Chef_lattice()
{
}

Chef_elements
Chef_lattice::get_chef_elements(Lattice_element & lattice_element)
{
	Chef_elements retval;
    Begin_end begin_end(element_map[&lattice_element]);
    for(int index = begin_end.begin; index < begin_end.end; ++index) {
    	retval.push_back(*beamline_iterators.at(index));
    }
    return retval;
}

Chef_lattice_section_sptr
Chef_lattice::get_chef_section_sptr(Chef_lattice_sptr this_chef_lattice_sptr,
        Lattice_element_slice & lattice_element_slice)
{
    if (this_chef_lattice_sptr.get() != this) {
        throw std::runtime_error(
                "get_chef_section_sptr requires a shared pointer to its parent instance in the first argument\n");
    }
    if (!have_sliced_beamline_) {
        throw std::runtime_error(
                "get_chef_section_sptr(Lattice_element_slice const&) called before construct_sliced_beamline\n");
    }
    if (element_slice_map.count(&lattice_element_slice) == 0) {
        throw std::runtime_error(
                "get_chef_section_sptr(Lattice_element_slice const&): slice not found\n");
    }
    Begin_end begin_end(element_slice_map[&lattice_element_slice]);
    return Chef_lattice_section_sptr(
            new Chef_lattice_section(this_chef_lattice_sptr, begin_end.begin,
                    begin_end.end));
}

ElmPtr
slice_chef_element(ElmPtr & elm, double left, double right, double tolerance)
{
    double length = elm->Length();
    ElmPtr retval, left_part, right_part;
    if (left == 0.0) {
        if (floating_point_equal(length, right, tolerance)) {
            retval = elm;
        } else {
            elm->Split(right / length, left_part, right_part);
            retval = left_part;
        }
    } else {
        elm->Split(left / length, left_part, right_part);
        if (floating_point_equal(length, right, tolerance)) {
            retval = right_part;
        } else {
            ElmPtr second_left_part, second_right_part;
            right_part->Split((right - left) / (length - left),
                    second_left_part, second_right_part);
            int index = int(left / (right - left)) + 1;
            std::stringstream element_name(stringstream::out);
            element_name << elm->Name() << "_1_" << index;
            second_left_part->rename(element_name.str().c_str());
            retval = second_left_part;
        }
    }

    return retval;
}

Chef_elements
Chef_lattice::get_chef_elements_from_slice(Lattice_element_slice const& slice)
{
    //std::cout << "\ngetting elements from slice: " << slice.as_string() << std::endl;
    Chef_elements all_elements = lattice_sptr->get_element_adaptor_map().get_adaptor(
            slice.get_lattice_element().get_type())->get_chef_elements(
            slice.get_lattice_element(), brho);
    Chef_elements retval;
    if (slice.is_whole()) {
        //std::cout << "slice is whole" << std::endl;
        retval = all_elements;
    } else {
        // the lattice element corresponds to a list of chef elements.  The lattice_element_slice
        // could correspond to a subset of them and even start midway through the list.  The first
        // order of business is to skip those chef elements that fully precede the slice.
        // std::cout << "slice is not whole, peeling off parts of elements" << std::endl;
        const double tolerance = 1.0e-8;
        double left = slice.get_left();
        double right = slice.get_right();
        double s = 0.0;
        //std::cout << "left: " << left << ", right: " << right << std::endl;
        Chef_elements::iterator c_it = all_elements.begin();
        bool complete = false;
        double element_left, element_right;
        double total_done = 0.0;
        while (!complete) {
            //std::cout << "not complete, considering current element: " << chef_element_as_string(*c_it) << std::endl;
            //std::cout << "total_done: " << std::setprecision(16) << total_done << ", need to reach " << right-left << std::endl;
            element_left = left - s + total_done;
            element_right = 0.0;
            // element_left current position relative to start of slice
            // s is the starting s position of this element relative to the start of the slice
            double chef_element_length = (*c_it)->Length();
            //std::cout << "s: " << std::setprecision(16) << s << ", element_left: "  << element_left << ", element_right: " << element_right << ", chef_element_length: " << chef_element_length << std::endl;
            if (floating_point_leq(chef_element_length, element_left, tolerance)
                    && (total_done == 0.0)) {
                //is this element fully contained before the slice?
                s += chef_element_length;
                // skip this chef element because the current slice starts after the element
                ++c_it;
                //std::cout << "need more length than this element" << std::endl;
                if (c_it == all_elements.end()) {
                    throw(std::runtime_error(
                            "get_chef_section_from_slice iterated beyond end of element list"));
                }
            } else {
                //std::cout << "looks like we need more than current element" << std::endl;
                if (chef_element_length == 0.0) {
                    retval.push_back(*c_it);
                    ++c_it;
                    //std::cout << "current element has 0 length, skipping" << std::endl;
                } else {
                    if (floating_point_leq(right, s + chef_element_length,
                            tolerance)) {
                        // take part of element
                        element_right = right - s;
                        //std::cout << "take slice of element [" << std::setprecision(16) << element_left << "," << element_right << "]" << std::endl;
                        retval.push_back(
                                slice_chef_element(*c_it, element_left,
                                        element_right, tolerance));
                        s += element_right - element_left;
                        total_done += element_right - element_left;
                        //std::cout << "total_done: " << std::setprecision(16) << total_done << std::endl;
                    } else {
                        // take rest of element
                        element_right = chef_element_length;
                        //std::cout << "take rest of element" << std::setprecision(16) << element_left << "," << element_right << "]" << std::endl;
                        retval.push_back(
                                slice_chef_element(*c_it, element_left,
                                        element_right, tolerance));
                        s += chef_element_length;
                        total_done += element_right - element_left;
                        //std::cout << "total_done: " << std::setprecision(16) << total_done << std::endl;
                        ++c_it;
                    }
                }
                if (floating_point_equal(total_done, right - left, tolerance)) {
                    complete = true;
                    if (floating_point_equal(element_right,
                            chef_element_length, tolerance)) {
                        while ((++c_it != all_elements.end())
                                && ((*c_it)->Length() == 0.0)) {
                            retval.push_back(*c_it);
                        }
                    }
                }
            }
        }
    }

    return retval;
}

Lattice_element &
Chef_lattice::get_lattice_element(ElmPtr const& chef_element)
{
    Lattice_element * lattice_element_ptr;
    bool found = false;
    for (std::map<Lattice_element*, Begin_end >::iterator it =
            element_map.begin(); it != element_map.end(); ++it) {
    	for (int index=it->second.begin; index != it->second.end; ++index) {
    		if ((*beamline_iterators.at(index)) == chef_element) {
                found = true;
                lattice_element_ptr = it->first;
                break;
            }
        }
        if (found) {
            break;
        }
    }
    if (!found) {
        throw std::runtime_error(
                "Chef_lattice::get_lattice_element: no match for chef element");
    }
    return *lattice_element_ptr;
}

Lattice_element_slice &
Chef_lattice::get_lattice_element_slice(ElmPtr const& chef_element)
{
    Lattice_element_slice * lattice_element_slice_ptr;
    bool found = false;
    for (std::map<Lattice_element_slice*, Begin_end >::iterator it =
            element_slice_map.begin(); it != element_slice_map.end(); ++it) {
    	for (int index=it->second.begin; index != it->second.end; ++index) {
    		if ((*sliced_beamline_iterators.at(index)) == chef_element) {
                found = true;
                lattice_element_slice_ptr = it->first;
                break;
            }
        }
        if (found) {
            break;
        }
    }
    if (!found) {
        throw std::runtime_error(
                "Chef_lattice::get_lattice_element_slice: no match for chef element");
    }
    return *lattice_element_slice_ptr;
}

bool
Chef_lattice::have_sliced_beamline() const
{
    return have_sliced_beamline_;
}

void
Chef_lattice::construct_sliced_beamline(Lattice_element_slices const& slices)
{
    this->slices = slices;
    BmlPtr unpolished_beamline_sptr(new beamline());
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        Chef_elements chef_elements = get_chef_elements_from_slice(*(*it));
        for (Chef_elements::const_iterator c_it = chef_elements.begin(); c_it
                != chef_elements.end(); ++c_it) {
            unpolished_beamline_sptr->append(*c_it);
        }
        unpolished_beamline_sptr->append(lattice_element_marker);
    }
    sliced_beamline_sptr = polish_beamline(unpolished_beamline_sptr);
    extract_element_slice_map(slices);
    lattice_element_slices = slices;
    have_sliced_beamline_ = true;
}

BmlPtr
Chef_lattice::get_beamline_sptr()
{
    return beamline_sptr;
}

BmlPtr
Chef_lattice::get_sliced_beamline_sptr()
{
    if (!have_sliced_beamline_) {
        throw std::runtime_error(
                "get_sliced_beamline_sptr() called before construct_sliced_beamline\n");
    }
    return sliced_beamline_sptr;
}

beamline::iterator
Chef_lattice::get_sliced_beamline_iterator(int index)
{
    return sliced_beamline_iterators.at(index);
}

beamline::const_iterator
Chef_lattice::get_sliced_beamline_const_iterator(int index) const
{
    return sliced_beamline_const_iterators.at(index);
}

template<class Archive>
    void
    Chef_lattice::save(Archive & ar, const unsigned int version) const
    {
        ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
        ar & BOOST_SERIALIZATION_NVP(slices);
        ar & BOOST_SERIALIZATION_NVP(have_sliced_beamline_);
        ar & BOOST_SERIALIZATION_NVP(brho);
    }

template<class Archive>
    void
    Chef_lattice::load(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
        ar & BOOST_SERIALIZATION_NVP(slices);
        ar & BOOST_SERIALIZATION_NVP(have_sliced_beamline_);
        ar & BOOST_SERIALIZATION_NVP(brho);
        lattice_element_marker = ElmPtr(
                new marker(internal_marker_name));
        construct_beamline();
        if (have_sliced_beamline_) {
            construct_sliced_beamline(slices);
        }
    }

template
void
Chef_lattice::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Chef_lattice::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Chef_lattice::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Chef_lattice::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Chef_lattice::~Chef_lattice()
{

}
