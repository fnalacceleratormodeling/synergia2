#ifndef CHEF_LATTICE_H_
#define CHEF_LATTICE_H_

#include <boost/enable_shared_from_this.hpp>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/chef_elements.h"
#include "synergia/lattice/chef_lattice_section_fwd.h"
#include <beamline/beamline_elements.h>

class Chef_lattice : public boost::enable_shared_from_this<Chef_lattice >
{
    friend class Chef_lattice_tester;
private:
    struct Begin_end
    {
        int begin, end;
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(begin);
                ar & BOOST_SERIALIZATION_NVP(end);
            }
    };

    Lattice_sptr lattice_sptr;
    Lattice_element_slices slices;
    Element_adaptor_map_sptr element_adaptor_map_sptr;
    BmlPtr beamline_sptr;
    BmlPtr sliced_beamline_sptr;
    Lattice_element_slices lattice_element_slices;
    bool have_sliced_beamline_;
    std::vector<beamline::iterator > sliced_beamline_iterators;
    std::vector<beamline::const_iterator > sliced_beamline_const_iterators;
    ElmPtr lattice_element_marker;
    std::map<const Lattice_element*, Begin_end > element_map;
    std::map<const Lattice_element_slice*, Begin_end > element_slice_map;
    double brho;

    void
    construct_beamline();
    void
    register_beamline(beamline & the_beamline);
    BmlPtr
    polish_beamline(BmlPtr beamline_sptr);
    void
    extract_element_map();
    void
    extract_element_slice_map(Lattice_element_slices const& slices);
    Chef_elements
    get_chef_elements_from_slice(Lattice_element_slice const& slice);
    void
    construct();
public:
    Chef_lattice(Lattice_sptr lattice_sptr);
    Chef_lattice(Lattice_sptr lattice_sptr,
            Element_adaptor_map_sptr element_adaptor_map_sptr);
    double
    get_brho() const;
    // Default constructor for serialization use only
    Chef_lattice();
    Element_adaptor_map_sptr
    get_element_adaptor_map_sptr();
    Chef_lattice_section_sptr
    get_chef_section_sptr(Lattice_element const& lattice_element);
    Chef_lattice_section_sptr
    get_chef_section_sptr(Lattice_element_slice const& lattice_element_slice);
    Lattice_element const&
    get_lattice_element(ElmPtr const& chef_element);
    Lattice_element_slice const&
    get_lattice_element_slice(ElmPtr const& chef_element);
    bool
    have_sliced_beamline() const;
    void
    construct_sliced_beamline(Lattice_element_slices const& slices);
    BmlPtr
    get_beamline_sptr();
    BmlPtr
    get_sliced_beamline_sptr();
    beamline::iterator
    get_sliced_beamline_iterator(int index);
    beamline::const_iterator
    get_sliced_beamline_const_iterator(int index) const;
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const
        {
            ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
            ar & BOOST_SERIALIZATION_NVP(slices);
            ar & BOOST_SERIALIZATION_NVP(element_adaptor_map_sptr);
            ar & BOOST_SERIALIZATION_NVP(have_sliced_beamline_);
            ar & BOOST_SERIALIZATION_NVP(brho);
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
            ar & BOOST_SERIALIZATION_NVP(slices);
            ar & BOOST_SERIALIZATION_NVP(element_adaptor_map_sptr);
            ar & BOOST_SERIALIZATION_NVP(have_sliced_beamline_);
            ar & BOOST_SERIALIZATION_NVP(brho);
            lattice_element_marker = ElmPtr(
                    new marker("synergia_lattice_element_marker"));
            construct_beamline();
            if (have_sliced_beamline_) {
                construct_sliced_beamline(slices);
            }
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    ~Chef_lattice();
    static const char internal_marker_name[];
};

typedef boost::shared_ptr<Chef_lattice > Chef_lattice_sptr;

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

#endif /* CHEF_LATTICE_H_ */
