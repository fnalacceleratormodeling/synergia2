#include "ff_drift.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"
#include <iomanip>

FF_drift::FF_drift()
{
}

double FF_drift::get_reference_cdt(double length,
                                   Reference_particle & reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(reference_particle.get_state()[Bunch::cdt]);
    double dpop(reference_particle.get_state()[Bunch::dpop]);
    double reference_momentum = reference_particle.get_momentum();
    double m = reference_particle.get_mass();

    double cdt_orig = cdt;
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m, 0.0);

    return cdt - cdt_orig;
}

void FF_drift::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    double length = slice.get_right() - slice.get_left();

    typedef PropagatorTraits<JetParticle>::State_t State_t;
    typedef PropagatorTraits<JetParticle>::Component_t Component_t;

    State_t& state = jet_particle.State();

    Component_t & x(state[Chef::x]);
    Component_t & xp(state[Chef::xp]);
    Component_t & y(state[Chef::y]);
    Component_t & yp(state[Chef::yp]);
    Component_t & cdt(state[Chef::cdt]);
    Component_t & dpop(state[Chef::dpop]);

    double reference_momentum = jet_particle.ReferenceMomentum();
    double m = jet_particle.Mass();

    Particle chef_particle(jet_particle);
    Reference_particle reference_particle(
                chef_particle_to_reference_particle(chef_particle));
    double reference_cdt = get_reference_cdt(length, reference_particle);

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
}

void FF_drift::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double t0 = MPI_Wtime();
     const double length = slice.get_right() - slice.get_left();
     const int local_num = bunch.get_local_num();
     const double reference_momentum = bunch.get_reference_particle().get_momentum();
     const double m = bunch.get_mass();
     const double reference_cdt = get_reference_cdt(length,
                                                    bunch.get_reference_particle());
     double * RESTRICT xa, * RESTRICT xpa, * RESTRICT ya, * RESTRICT ypa,
             * RESTRICT cdta, * RESTRICT dpopa;
     bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

     double xtmp[local_num];
     double ytmp[local_num];
     double cdttmp[local_num];

     const int num_blocks = local_num / GSVector::size;
     const int block_last = num_blocks * GSVector::size;
     double t1 = MPI_Wtime();
     #pragma omp parallel for
     for (int part = 0; part < block_last; part += GSVector::size) {
         GSVector x(&xa[part]);
         GSVector xp(&xpa[part]);
         GSVector y(&ya[part]);
         GSVector yp(&ypa[part]);
         GSVector cdt(&cdta[part]);
         GSVector dpop(&dpopa[part]);

         FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
                    reference_cdt);

//         x.store(&xa[part]);
//         y.store(&ya[part]);
//         cdt.store(&cdta[part]);
         x.store(&xtmp[part]);
         y.store(&ytmp[part]);
         cdt.store(&cdttmp[part]);
     }
     #pragma omp parallel for
     for (int part = block_last; part < local_num; ++part) {
         double x(xa[part]);
         double xp(xpa[part]);
         double y(ya[part]);
         double yp(ypa[part]);
         double cdt(cdta[part]);
         double dpop(dpopa[part]);

         FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
                    reference_cdt);

         xa[part] = x;
         ya[part] = y;
         cdta[part] = cdt;
     }
     double t2 = MPI_Wtime();
//    #pragma omp parallel for
//     for(int part = 0; part < block_last; ++part) {
//         xa[part] = xtmp[part];
//         ya[part] = ytmp[part];
//         cdta[part] = cdttmp[part];
//     }
     std::copy(xtmp, xtmp+block_last, xa);
     std::copy(ytmp, ytmp+block_last, ya);
     std::copy(cdttmp, cdttmp+block_last, cdta);
     double t3 = MPI_Wtime();
     Logger logger(0);
//     logger << "jfa: GSVector::implentation " << GSVector::implementation << std::endl;
     logger << std::setw(8) << std::setprecision(6);
     logger << "drift-time: " << t1 -t0 << ", " << t2-t1 << ", " << t3-t2 << std::endl;
    bunch.get_reference_particle().increment_trajectory(length);
}

template<class Archive>
    void
    FF_drift::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_drift::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_drift::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_drift::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_drift::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_drift::~FF_drift()
{

}
