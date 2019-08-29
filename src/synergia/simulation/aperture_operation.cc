#include "synergia/simulation/aperture_operation.h"


#if 0
namespace aperture_impl
{
    template<class AP>
    struct discard_checker
    {
        typedef int value_type;

        AP ap;
        ConstParticles parts;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& discarded) const
        {
            if (ap.discard(parts, i))
            {
                //discard[i] = 1;
                ++discarded;
            }
            else
            {
                //discard[i] = 0;
            }
        }
    };

    struct particle_mover
    {
        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
        }
    };
}


void 
Aperture_operation::apply_impl(Bunch & bunch, Logger & logger) const
{
    using namespace aperture_impl;

    int discarded = 0;
    discard_checker<AP> dc;
    Kokkos::parallel_reduce(nparts, aa, discarded);

    particle_mover pm;
    Kokkos::parallel(1, pm);
}
#endif




#if 0
    double t0 = MPI_Wtime();
    bool write_loss=false;
    int b_index=-1; // AM: this value is written in the aperture_loss file when the bunch has no bucket index assigned
    Diagnostics_losses diagnostics_list=
         get_slice_sptr()->get_lattice_element().get_lattice().get_loss_diagnostics_list();
    Diagnostics_loss_sptr diagnostics_sptr;
    for (Diagnostics_losses::const_iterator d_it = diagnostics_list.begin();
        d_it != diagnostics_list.end(); ++d_it){
            if ( ((*d_it)->get_bunch().get_bucket_index()==bunch.get_bucket_index()) &&
                ((*d_it)->get_type()==Diagnostics_loss::aperture_type) )
            { 
              diagnostics_sptr=(*d_it);
              write_loss=true;
            }
    }          

    int nt;
    #pragma omp parallel
    { nt = omp_get_num_threads(); }

    if   (bunch.is_bucket_index_assigned())  b_index=bunch.get_bucket_index();
    int repetition=bunch.get_reference_particle().get_repetition();
    double s=bunch.get_reference_particle().get_s();
    double s_n=bunch.get_reference_particle().get_s_n();
    MArray1d coords(boost::extents[7]);

    MArray2d_ref particles(bunch.get_local_particles());
    MArray2d_ref s_particles(bunch.get_local_spectator_particles());

    int npart = bunch.get_local_num();
    int npart_s = bunch.get_local_spectator_num();

    int * discard = new int[npart];
    int * discard_count = new int[nt];

    int * discard_s = new int[npart_s];
    int * discard_s_count = new int[nt];

    int part_per_thread = npart / nt;
    int s_part_per_thread = npart_s / nt;

    #pragma omp parallel shared(nt, npart, npart_s, particles, s_particles, discard, discard_s, discard_count, discard_s_count)
    {
        int it = omp_get_thread_num();

        discard_count[it] = 0;
        discard_s_count[it] = 0;

        int s = it * part_per_thread;
        int e = (it==nt-1) ? npart : (s+part_per_thread);

        for (int part = s; part < e; ++part)
        {
            if (t(particles, part)) 
            {
                discard[part] = 1;
                ++discard_count[it];
            }
            else
            {
                discard[part] = 0;
            }
        }

        s = it * s_part_per_thread;
        e = (it==nt-1) ? npart_s : (s + s_part_per_thread);

        for (int part = s; part < e; ++part)
        {
            if (t(s_particles, part)) 
            {
                discard_s[part] = 1;
                ++discard_s_count[it];
            }
            else
            {
                discard_s[part] = 0;
            }
        }
    }

    // number of discarded particles
    int discarded = 0;
    int discarded_s = 0;

    for (int i=0; i<nt; ++i) 
    {
        discarded += discard_count[i];
        discarded_s += discard_s_count[i];
    }

    // arrange the particle array
    {
        // move all the discarded particles to the tail
        int head = 0;
        int tail = npart - 1;

        do
        {
            while (!discard[head] && head<tail) ++head;
            if (head >= tail) break;

            while ( discard[tail] && tail>head) --tail;
            if (head >= tail) break;

            double p0 = particles[head][0];
            double p1 = particles[head][1];
            double p2 = particles[head][2];
            double p3 = particles[head][3];
            double p4 = particles[head][4];
            double p5 = particles[head][5];
            double p6 = particles[head][6];

            particles[head][0] = particles[tail][0];
            particles[head][1] = particles[tail][1];
            particles[head][2] = particles[tail][2];
            particles[head][3] = particles[tail][3];
            particles[head][4] = particles[tail][4];
            particles[head][5] = particles[tail][5];
            particles[head][6] = particles[tail][6];

            particles[tail][0] = p0;
            particles[tail][1] = p1;
            particles[tail][2] = p2;
            particles[tail][3] = p3;
            particles[tail][4] = p4;
            particles[tail][5] = p5;
            particles[tail][6] = p6;

            ++head;
            --tail;

        } while(head < tail);

        // move some lost particles over to the padding area
        int padded  = bunch.get_local_num_padded();
        int padding = padded - npart;
        int np = discarded < padding ? discarded : padding;

        for (int i=0; i<np; ++i)
        {
            // pl: position of next lost particle
            // pp: position of next padding slot
            int pl = npart - discarded + i;
            int pp = padded - 1 - i;

            // copy the lost particle over to the padding slot
            particles[pp][0] = particles[pl][0];
            particles[pp][1] = particles[pl][1];
            particles[pp][2] = particles[pl][2];
            particles[pp][3] = particles[pl][3];
            particles[pp][4] = particles[pl][4];
            particles[pp][5] = particles[pl][5];
            particles[pp][6] = particles[pl][6];

            // makes pl the new padding slot
            particles[pl][0] = 0.0;
            particles[pl][1] = 0.0;
            particles[pl][2] = 0.0;
            particles[pl][3] = 0.0;
            particles[pl][4] = 0.0;
            particles[pl][5] = 0.0;
            particles[pl][6] = 0.0;
        }

        // finalize the bunch for new particle array pointers
        double charge = (discarded > 0) ? discarded * bunch.get_real_num() / bunch.get_total_num() : 0.0;
        deposit_charge(charge);
        bunch.set_local_num(npart - discarded);
    }

    // arrange the spectator particle array
    {
        // move all the discarded spectator particles to the tail
        int head = 0;
        int tail = npart_s - 1;

        do
        {
            while (!discard_s[head] && head<tail) ++head;
            if (head >= tail) break;

            while ( discard_s[tail] && tail>head) --tail;
            if (head >= tail) break;

            double p0 = s_particles[head][0];
            double p1 = s_particles[head][1];
            double p2 = s_particles[head][2];
            double p3 = s_particles[head][3];
            double p4 = s_particles[head][4];
            double p5 = s_particles[head][5];
            double p6 = s_particles[head][6];

            s_particles[head][0] = s_particles[tail][0];
            s_particles[head][1] = s_particles[tail][1];
            s_particles[head][2] = s_particles[tail][2];
            s_particles[head][3] = s_particles[tail][3];
            s_particles[head][4] = s_particles[tail][4];
            s_particles[head][5] = s_particles[tail][5];
            s_particles[head][6] = s_particles[tail][6];

            s_particles[tail][0] = p0;
            s_particles[tail][1] = p1;
            s_particles[tail][2] = p2;
            s_particles[tail][3] = p3;
            s_particles[tail][4] = p4;
            s_particles[tail][5] = p5;
            s_particles[tail][6] = p6;

            ++head;
            --tail;

        } while(head < tail);

        // move some lost spectator particles over to the padding area
        int padded  = bunch.get_local_spectator_num_padded();
        int padding = padded - npart_s;
        int np = discarded_s < padding ? discarded_s : padding;

        for (int i=0; i<np; ++i)
        {
            // pl: position of next lost particle
            // pp: position of next padding slot
            int pl = npart - discarded + i;
            int pp = padded - 1 - i;

            // copy the lost particle over to the padding slot
            s_particles[pp][0] = s_particles[pl][0];
            s_particles[pp][1] = s_particles[pl][1];
            s_particles[pp][2] = s_particles[pl][2];
            s_particles[pp][3] = s_particles[pl][3];
            s_particles[pp][4] = s_particles[pl][4];
            s_particles[pp][5] = s_particles[pl][5];
            s_particles[pp][6] = s_particles[pl][6];

            // makes pl the new padding slot
            s_particles[pl][0] = 0.0;
            s_particles[pl][1] = 0.0;
            s_particles[pl][2] = 0.0;
            s_particles[pl][3] = 0.0;
            s_particles[pl][4] = 0.0;
            s_particles[pl][5] = 0.0;
            s_particles[pl][6] = 0.0;
        }

        bunch.set_local_spectator_num(npart_s - discarded_s);
    }

    double t1 = MPI_Wtime();
    if (verbosity > 5) 
    {
        logger << "Aperture_operation: type = " << get_aperture_type()
               << ", discarded: " << discarded
               << ", discarded spectators: " << discarded_s
               << ", time = " << std::fixed << std::setprecision(3) << t1
                - t0 << "s_n" << std::endl;
    }

    delete [] discard;
    delete [] discard_count;

    delete [] discard_s;
    delete [] discard_s_count;
}
#endif


#if 0
const double Circular_aperture_operation::default_radius(1000.0);
const char Circular_aperture_operation::aperture_type[] = "circular";
const char Circular_aperture_operation::attribute_name[] = "circular";

Circular_aperture_operation::Circular_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{

    if (slice_sptr->get_lattice_element().has_double_attribute(
            "circular_aperture_radius")) {
        radius = slice_sptr->get_lattice_element().get_double_attribute(
                "circular_aperture_radius");
    } else {
        radius = default_radius;
    }
    radius2 = radius * radius;
}

Circular_aperture_operation::Circular_aperture_operation()
{
}

const char *
Circular_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Circular_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Circular_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Circular_aperture_operation::operator==(
        Circular_aperture_operation const& circular_aperture_operation) const
{
    return (radius == circular_aperture_operation.radius);
}

void
Circular_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "circular_aperture-apply");
}

template<class Archive>
    void
    Circular_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(radius);
        ar & BOOST_SERIALIZATION_NVP(radius2);
    }

template
void
Circular_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Circular_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Circular_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Circular_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Circular_aperture_operation::~Circular_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Circular_aperture_operation);

const char Elliptical_aperture_operation::aperture_type[] = "elliptical";
const char Elliptical_aperture_operation::attribute_name[] = "elliptical";

Elliptical_aperture_operation::Elliptical_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "elliptical_aperture_horizontal_radius")) {
        horizontal_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_horizontal_radius attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "elliptical_aperture_vertical_radius")) {
        vertical_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Elliptical_aperture_operation: elliptical_aperture requires an elliptical_aperture_vertical_radius attribute");
    }
    h2 = horizontal_radius * horizontal_radius;
    v2 = vertical_radius * vertical_radius;
}

Elliptical_aperture_operation::Elliptical_aperture_operation()
{
}

const char *
Elliptical_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Elliptical_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Elliptical_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Elliptical_aperture_operation::operator==(
        Elliptical_aperture_operation const& elliptical_aperture_operation) const
{
    return ((horizontal_radius
            == elliptical_aperture_operation.horizontal_radius)
            && (vertical_radius
                    == elliptical_aperture_operation.vertical_radius));
}

void
Elliptical_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "elliptical_aperture-apply");
}

template<class Archive>
    void
    Elliptical_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(horizontal_radius);
        ar & BOOST_SERIALIZATION_NVP(vertical_radius);
        ar & BOOST_SERIALIZATION_NVP(h2);
        ar & BOOST_SERIALIZATION_NVP(v2);
    }

template
void
Elliptical_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Elliptical_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Elliptical_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Elliptical_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Elliptical_aperture_operation::~Elliptical_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Elliptical_aperture_operation);

const char Rectangular_aperture_operation::aperture_type[] = "rectangular";
const char Rectangular_aperture_operation::attribute_name[] = "rectangular";

Rectangular_aperture_operation::Rectangular_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_width")) {
        width = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_width");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_width attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_height")) {
        height = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_height");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_height attribute");
    }
}

Rectangular_aperture_operation::Rectangular_aperture_operation()
{
}

const char *
Rectangular_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Rectangular_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Rectangular_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Rectangular_aperture_operation::operator==(
        Rectangular_aperture_operation const& rectangular_aperture_operation) const
{
    return ((width == rectangular_aperture_operation.width) && (height
            == rectangular_aperture_operation.height));
}

void
Rectangular_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "rectangular_aperture-apply");
}

template<class Archive>
    void
    Rectangular_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(width);
        ar & BOOST_SERIALIZATION_NVP(height);
    }

template
void
Rectangular_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rectangular_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rectangular_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rectangular_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rectangular_aperture_operation::~Rectangular_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rectangular_aperture_operation);

const char Rectangular_with_ears_aperture_operation::aperture_type[] = "rectangular_with_ears";
const char Rectangular_with_ears_aperture_operation::attribute_name[] = "rectangular_with_ears";

Rectangular_with_ears_aperture_operation::Rectangular_with_ears_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_width")) {
        width = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_width");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_width attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_height")) {
        height = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_height");
    } else {
        throw std::runtime_error(
                "Rectangular_aperture_operation: rectangular_aperture requires an rectangular_aperture_height attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "rectangular_aperture_ear_offset")) {
        ear_offset = slice_sptr->get_lattice_element().get_double_attribute(
                "rectangular_aperture_ear_offset");
    } else {
        throw std::runtime_error(
                "Rectangular_with_ears_aperture_operation: rectangular_with_ears_aperture requires an rectangular_aperture_ear_offset attribute");
    }
    radius = 0.5*height - ear_offset;
}

Rectangular_with_ears_aperture_operation::Rectangular_with_ears_aperture_operation()
{
}

const char *
Rectangular_with_ears_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Rectangular_with_ears_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Rectangular_with_ears_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Rectangular_with_ears_aperture_operation::operator==(
        Rectangular_with_ears_aperture_operation const& rectangular_with_ears_aperture_operation) const
{
    return ((width == rectangular_with_ears_aperture_operation.width) &&
            (height == rectangular_with_ears_aperture_operation.height) &&
            (ear_offset == rectangular_with_ears_aperture_operation.ear_offset));
}

void
Rectangular_with_ears_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "rectangular_aperture-apply");
}

template<class Archive>
    void
    Rectangular_with_ears_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(width);
        ar & BOOST_SERIALIZATION_NVP(height);
        ar & BOOST_SERIALIZATION_NVP(ear_offset);
        ar & BOOST_SERIALIZATION_NVP(radius);
    }

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Rectangular_with_ears_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Rectangular_with_ears_aperture_operation::~Rectangular_with_ears_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Rectangular_with_ears_aperture_operation);

const char Polygon_aperture_operation::aperture_type[] = "polygon";
const char Polygon_aperture_operation::attribute_name[] = "polygon";

Polygon_aperture_operation::Polygon_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "the_number_of_vertices")) {
        num_vertices = int(
                slice_sptr->get_lattice_element().get_double_attribute(
                        "the_number_of_vertices"));
        if (num_vertices < 3) throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires at least 3 vertices");
    } else {
        throw std::runtime_error(
                "Polygon_aperture_operation: polygon_aperture requires the_number_of vertices attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "min_radius2")) {
        min_radius2 = slice_sptr->get_lattice_element().get_double_attribute(
                        "min_radius2");
    } else {
        min_radius2 = 0.0;
    }
    for (int index = 0; index < num_vertices; ++index) {
        std::string ss = boost::lexical_cast<std::string >(index + 1);
        std::string x = "pax" + ss;
        std::string y = "pay" + ss;
        if ((slice_sptr->get_lattice_element().has_double_attribute(x))
                && (slice_sptr->get_lattice_element().has_double_attribute(y))) {
            vertices.push_back(
                    std::complex<double >(
                            slice_sptr->get_lattice_element().get_double_attribute(
                                    x),
                            slice_sptr->get_lattice_element().get_double_attribute(
                                    y)));
        } else {
            throw std::runtime_error(
                    "Polygon_aperture_operation: polygon_aperture requires x and y coordinate attributes for each vertex");
        }
    }
}

Polygon_aperture_operation::Polygon_aperture_operation()
{
}

const char *
Polygon_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Polygon_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Polygon_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Polygon_aperture_operation::operator==(
        Polygon_aperture_operation const& polygon_aperture_operation) const
{
    return ((vertices == polygon_aperture_operation.vertices)
            && (min_radius2 == polygon_aperture_operation.min_radius2));
}

void
Polygon_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "polygon_aperture-apply");
}

template<class Archive>
    void
    Polygon_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(num_vertices);
        ar & BOOST_SERIALIZATION_NVP(vertices);
        ar & BOOST_SERIALIZATION_NVP(min_radius2);
    }

template
void
Polygon_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Polygon_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Polygon_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Polygon_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Polygon_aperture_operation::~Polygon_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Polygon_aperture_operation);

const char Wire_elliptical_aperture_operation::aperture_type[] =
        "wire_elliptical";
const char Wire_elliptical_aperture_operation::attribute_name[] =
        "wire_elliptical";

Wire_elliptical_aperture_operation::Wire_elliptical_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_horizontal_radius")) {
        horizontal_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "wire_elliptical_aperture_horizontal_radius");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_horizontal_radius attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_vertical_radius")) {
        vertical_radius
                = slice_sptr->get_lattice_element().get_double_attribute(
                        "wire_elliptical_aperture_vertical_radius");
    } else {
        throw std::runtime_error(
                "Wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_vertical_radius attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_wire_x")) {
        wire_x = slice_sptr->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_wire_x");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_x attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_wire_width")) {
        wire_width = slice_sptr->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_wire_width");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_wire_width attribute");
    }
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "wire_elliptical_aperture_gap")) {
        gap = slice_sptr->get_lattice_element().get_double_attribute(
                "wire_elliptical_aperture_gap");
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_elliptical_aperture requires an wire_elliptical_aperture_gap attribute");
    }
    h2 = horizontal_radius * horizontal_radius;
    v2 = vertical_radius * vertical_radius;
}

Wire_elliptical_aperture_operation::Wire_elliptical_aperture_operation()
{
}

const char *
Wire_elliptical_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Wire_elliptical_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Wire_elliptical_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Wire_elliptical_aperture_operation::operator==(
        Wire_elliptical_aperture_operation const& wire_elliptical_aperture_operation) const
{
    return ((horizontal_radius
            == wire_elliptical_aperture_operation.horizontal_radius)
            && (vertical_radius
                    == wire_elliptical_aperture_operation.vertical_radius)
            && (wire_x == wire_elliptical_aperture_operation.wire_x)
            && (wire_width == wire_elliptical_aperture_operation.wire_width)
            && (gap == wire_elliptical_aperture_operation.gap));
}

void
Wire_elliptical_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "wire_elliptical_aperture-apply");
}

template<class Archive>
    void
    Wire_elliptical_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(horizontal_radius);
        ar & BOOST_SERIALIZATION_NVP(vertical_radius);
        ar & BOOST_SERIALIZATION_NVP(h2);
        ar & BOOST_SERIALIZATION_NVP(v2);
        ar & BOOST_SERIALIZATION_NVP(wire_x);
        ar & BOOST_SERIALIZATION_NVP(wire_width);
        ar & BOOST_SERIALIZATION_NVP(gap);
    }

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Wire_elliptical_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Wire_elliptical_aperture_operation::~Wire_elliptical_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Wire_elliptical_aperture_operation);

const char Lambertson_aperture_operation::aperture_type[] = "lambertson";
const char Lambertson_aperture_operation::attribute_name[] = "lambertson";

Lambertson_aperture_operation::Lambertson_aperture_operation(
        Lattice_element_slice_sptr slice_sptr) :
    Aperture_operation(slice_sptr)
{
    if (slice_sptr->get_lattice_element().has_double_attribute(
            "lambertson_aperture_radius")) {
        radius = slice_sptr->get_lattice_element().get_double_attribute(
                "lambertson_aperture_radius");
    } else {
        throw std::runtime_error(
                "lambertson_aperture_operation: lambertson_aperture requires an lambertson_aperture_radius attribute");
    }
}

Lambertson_aperture_operation::Lambertson_aperture_operation()
{
}

const char *
Lambertson_aperture_operation::get_aperture_type() const
{
    return aperture_type;
}

bool
Lambertson_aperture_operation::operator==(
        Aperture_operation const& aperture_operation) const
{
    if (aperture_type == aperture_operation.get_aperture_type()) {
        return operator==(
                *static_cast<Lambertson_aperture_operation const* > (&aperture_operation));
    } else {
        return false;
    }
}

bool
Lambertson_aperture_operation::operator==(
        Lambertson_aperture_operation const& lambertson_aperture_operation) const
{
    return (radius == lambertson_aperture_operation.radius);
}

void
Lambertson_aperture_operation::apply(Bunch & bunch, int verbosity, Logger & logger)
{
    double t;
    //t = simple_timer_current();
    dump_particles(*this, bunch, verbosity, logger);
    //apply_impl(*this, bunch, verbosity, logger);
    //t = simple_timer_show(t, "lambertson_aperture-apply");
}

template<class Archive>
    void
    Lambertson_aperture_operation::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Aperture_operation);
        ar & BOOST_SERIALIZATION_NVP(radius);
    }

template
void
Lambertson_aperture_operation::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lambertson_aperture_operation::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);
    
template
void
Lambertson_aperture_operation::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lambertson_aperture_operation::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Lambertson_aperture_operation::~Lambertson_aperture_operation()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Lambertson_aperture_operation)
#endif
