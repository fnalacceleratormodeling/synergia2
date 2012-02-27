#ifndef LATTICE_SIMULATOR_H_
#define LATTICE_SIMULATOR_H_

#include <string>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/simulation/operation_extractor.h"
#include "synergia/simulation/aperture_operation_extractor.h"
#include "synergia/simulation/step.h"
#include <physics_toolkit/LattFuncSage.h>
#include <physics_toolkit/BeamlineContext.h>
#include <physics_toolkit/normalFormSage.h>

#include <string>
#include "synergia/utils/serialization.h"

typedef boost::shared_ptr<normalFormSage> Normal_form_sage_sptr;

struct Lattice_functions
{
    Lattice_functions();
    Lattice_functions(LattFuncSage::lattFunc const& latt_func);
    double alpha_x, alpha_y;
    double beta_x, beta_y;
    double psi_x, psi_y;
    double D_x, D_y;
    double Dprime_x, Dprime_y;
    double arc_length;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(alpha_x);
            ar & BOOST_SERIALIZATION_NVP(alpha_y);
            ar & BOOST_SERIALIZATION_NVP(beta_x);
            ar & BOOST_SERIALIZATION_NVP(beta_y);
            ar & BOOST_SERIALIZATION_NVP(psi_x);
            ar & BOOST_SERIALIZATION_NVP(psi_y);
            ar & BOOST_SERIALIZATION_NVP(D_x);
            ar & BOOST_SERIALIZATION_NVP(D_y);
            ar & BOOST_SERIALIZATION_NVP(Dprime_x);
            ar & BOOST_SERIALIZATION_NVP(Dprime_y);
            ar & BOOST_SERIALIZATION_NVP(arc_length);
        }
};

class Lattice_simulator
{
private:
    Lattice_sptr lattice_sptr;
    Lattice_element_slices slices;
    bool have_slices;
    Chef_lattice_sptr chef_lattice_sptr;
    Operation_extractor_map_sptr extractor_map_sptr;
    Aperture_operation_extractor_map_sptr aperture_extractor_map_sptr;
    bool have_beamline_context;
    BmlContextPtr beamline_context_sptr;
    int map_order;
    double bucket_length;
    bool have_element_lattice_functions;
    bool have_slice_lattice_functions;
    double horizontal_tune, vertical_tune;
    bool have_tunes;
    std::map<Lattice_element *, Lattice_functions >
            lattice_functions_element_map;
    std::map<Lattice_element_slice *, Lattice_functions >
            lattice_functions_slice_map;
    MArray2d linear_one_turn_map;
	void
	construct_extractor_map();
    void
    construct_aperture_extractor_map();
    void
    calculate_beamline_context();
    BmlContextPtr
    get_beamline_context();
    void
    construct_sliced_chef_beamline();
    void
    get_tunes();
    Normal_form_sage_sptr normal_form_sage_sptr;
public:
    Lattice_simulator(Lattice_sptr lattice, int map_order);
    // Default constructor for serialization use only
    Lattice_simulator();
    void
    set_slices(Lattice_element_slices const& slices);
    int
    get_map_order() const;
    void
    set_bucket_length();
    /// bucket length is in z_lab frame
    double
    get_bucket_length();
    int
    get_number_buckets();
    Operation_extractor_map_sptr
    get_operation_extractor_map_sptr();
    Aperture_operation_extractor_map_sptr
    get_aperture_operation_extractor_map_sptr();
    Lattice_sptr
    get_lattice_sptr();
    Chef_lattice_sptr
    get_chef_lattice_sptr();
    void
    update();
    void
    calculate_element_lattice_functions();
    void
    calculate_slice_lattice_functions();
    Lattice_functions const&
    get_lattice_functions(Lattice_element & lattice_element);
    Lattice_functions const&
    get_lattice_functions(Lattice_element_slice & lattice_element_slice);
    double
    get_horizontal_tune();
    double
    get_vertical_tune();
    void calculate_normal_form();
    Normal_form_sage_sptr get_normal_form_sptr();
    Const_MArray2d_ref get_linear_one_turn_map();
    void convert_human_to_normal(MArray2d_ref coords);
    void convert_normal_to_human(MArray2d_ref coords);
    bool check_linear_normal_form();
    std::vector<double> get_stationary_actions(const double stdx, const double stdy, const double stdz);
    void
    adjust_tunes(double horizontal_tune, double vertical_tune,
            Lattice_elements const& horizontal_correctors,
            Lattice_elements const& vertical_correctors,
            double tolerance = 1.0e-6);
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const
        {
            ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
            ar & BOOST_SERIALIZATION_NVP(slices);
            ar & BOOST_SERIALIZATION_NVP(have_slices);
            ar & BOOST_SERIALIZATION_NVP(chef_lattice_sptr);
            ar & BOOST_SERIALIZATION_NVP(extractor_map_sptr);
            ar & BOOST_SERIALIZATION_NVP(aperture_extractor_map_sptr);
            ar & BOOST_SERIALIZATION_NVP(have_beamline_context);
            ar & BOOST_SERIALIZATION_NVP(map_order);
            ar & BOOST_SERIALIZATION_NVP(bucket_length);
            ar & BOOST_SERIALIZATION_NVP(have_element_lattice_functions);
            ar & BOOST_SERIALIZATION_NVP(have_slice_lattice_functions);
            ar & BOOST_SERIALIZATION_NVP(horizontal_tune);
            ar & BOOST_SERIALIZATION_NVP(vertical_tune);
            ar & BOOST_SERIALIZATION_NVP(have_tunes);
            ar & BOOST_SERIALIZATION_NVP(lattice_functions_element_map);
            ar & BOOST_SERIALIZATION_NVP(lattice_functions_slice_map);
            ar & BOOST_SERIALIZATION_NVP(linear_one_turn_map);
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
            ar & BOOST_SERIALIZATION_NVP(slices);
            ar & BOOST_SERIALIZATION_NVP(have_slices);
            ar & BOOST_SERIALIZATION_NVP(chef_lattice_sptr);
            ar & BOOST_SERIALIZATION_NVP(extractor_map_sptr);
            ar & BOOST_SERIALIZATION_NVP(aperture_extractor_map_sptr);
            ar & BOOST_SERIALIZATION_NVP(have_beamline_context);
            ar & BOOST_SERIALIZATION_NVP(map_order);
            ar & BOOST_SERIALIZATION_NVP(bucket_length);
            ar & BOOST_SERIALIZATION_NVP(have_element_lattice_functions);
            ar & BOOST_SERIALIZATION_NVP(have_slice_lattice_functions);
            ar & BOOST_SERIALIZATION_NVP(horizontal_tune);
            ar & BOOST_SERIALIZATION_NVP(vertical_tune);
            ar & BOOST_SERIALIZATION_NVP(have_tunes);
            ar & BOOST_SERIALIZATION_NVP(lattice_functions_element_map);
            ar & BOOST_SERIALIZATION_NVP(lattice_functions_slice_map);
            ar & BOOST_SERIALIZATION_NVP(linear_one_turn_map);
            if (have_beamline_context) {
                calculate_beamline_context();
            }
            normal_form_sage_sptr.reset();
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    ~Lattice_simulator();
};

#endif /* LATTICE_SIMULATOR_H_ */
