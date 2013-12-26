#include "bunch.h"
#include "synergia/utils/parallel_utils.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

const int Bunch::x;
const int Bunch::xp;
const int Bunch::y;
const int Bunch::yp;
const int Bunch::z;
const int Bunch::zp;
const int Bunch::cdt;
const int Bunch::dpop;
const int Bunch::id;


class Particle_id_offset
{
private:
    int offset;
public:
    Particle_id_offset() :
        offset(0)
    {
    }

    int
    get(int request_num, Commxx const & comm)
    {
        MPI_Bcast((void *) &offset, 1, MPI_INT, 0, comm.get());
        int old_offset = offset;
        int total_num;
        MPI_Reduce((void*) &request_num, (void*) &total_num, 1, MPI_INT,
                MPI_SUM, 0, comm.get());
        offset += total_num;
        return old_offset;
    }

};


static Particle_id_offset particle_id_offset;

void
Bunch::assign_ids(int local_offset)
{
    int global_offset, request_num;
    if (comm_sptr->get_rank() == 0) {
        request_num = total_num;
    } else {
        request_num = 0;
    }
    global_offset = particle_id_offset.get(request_num, *comm_sptr);
    for (int i = 0; i < local_num; ++i) {
        (*local_particles)[i][id] = i + local_offset + global_offset;
    }
}



template<class T, size_t C, int I>
    struct Sortable2d
    {
        typedef T data_type;
    };

template<class T, size_t C, int I>
    struct Sortable2d<T*, C, I >
    {
        typedef T data_type;
        typedef T arr_type[][C];
        typedef T row_type[C];
        struct Row
        {
            row_type data;
        };
        typedef Row cols_type[];

        Sortable2d(double* t, size_t sz) :
            ptr_(t), rows_(sz)
        {
        }

        struct Less
        {
            bool
            operator()(Row const& a, Row const& b)
            {
                return a.data[I] < b.data[I];
            }
        };

        Row*
        begin()
        {
            return (Row*) ptr_;
        }

        Row*
        end()
        {
            return (Row*) (ptr_ + (rows_ * C));
        }

        double* ptr_;
        size_t rows_;
    };

std::string
Bunch::get_local_particles_serialization_path() const
{
    std::stringstream sstream;
    sstream << "local_particles_";
    sstream << bucket_index;
    sstream << ".h5";
    return get_serialization_path(sstream.str());
}

void
Bunch::construct(int particle_charge, int total_num, double real_num)
{
    sort_counter = 0;
    sort_period = 10000;
    this->particle_charge = particle_charge;
    this->total_num = total_num;
    this->real_num = real_num;
    state = fixed_z;
    converter_ptr = &default_converter;
    if (comm_sptr->has_this_rank()) {
        local_num = decompose_1d_local(*comm_sptr, total_num);
        std::vector<int > offsets(comm_sptr->get_size()), counts(
                comm_sptr->get_size());
        decompose_1d(*comm_sptr, total_num, offsets, counts);
        local_num = counts[comm_sptr->get_rank()];
        local_particles = new MArray2d(boost::extents[local_num][7]);
        assign_ids(offsets[comm_sptr->get_rank()]);
    } else {
        local_num = 0;
        local_particles = new MArray2d(boost::extents[local_num][7]);
    }
}

Bunch::Bunch(Reference_particle const& reference_particle, int total_num,
        double real_num, Commxx_sptr comm_sptr) :
        z_period_length(0.0), z_periodic(0), reference_particle(
                reference_particle), bucket_index(0), comm_sptr(comm_sptr), default_converter()
{
    construct(reference_particle.get_charge(), total_num, real_num);
}

Bunch::Bunch(Reference_particle const& reference_particle, int total_num,
        double real_num, Commxx_sptr comm_sptr, int particle_charge) :
        z_period_length(0.0), z_periodic(0), reference_particle(
                reference_particle), bucket_index(0), comm_sptr(comm_sptr), default_converter()
{
    construct(particle_charge, total_num, real_num);
}

Bunch::Bunch(Reference_particle const& reference_particle, int total_num,
        double real_num, Commxx_sptr comm_sptr, double z_period_length,
        int bucket_index) :
        z_period_length(z_period_length), z_periodic(true), reference_particle(
                reference_particle), bucket_index(bucket_index), comm_sptr(
                comm_sptr), default_converter()
{
    construct(reference_particle.get_charge(), total_num, real_num);
}

Bunch::Bunch()
{
}


Bunch::Bunch(Bunch const& bunch) :
    reference_particle(bunch.reference_particle), comm_sptr(bunch.comm_sptr),
            default_converter()
{
    particle_charge = bunch.particle_charge;
    total_num = bunch.total_num;
    real_num = bunch.real_num;
    local_num = bunch.local_num;
    bucket_index=bunch.bucket_index;
    local_particles = new MArray2d(*(bunch.local_particles));
    state = bunch.state;
    z_period_length=bunch.z_period_length;
    z_periodic=bunch.z_periodic;
    if (bunch.converter_ptr == &(bunch.default_converter)) {
        converter_ptr = &default_converter;
    } else {
        converter_ptr = bunch.converter_ptr;
    }
}

Bunch &
Bunch::operator=(Bunch const& bunch)
{
    if (this != &bunch) {
        reference_particle = bunch.reference_particle;
        comm_sptr = bunch.comm_sptr;
        particle_charge = bunch.particle_charge;
        total_num = bunch.total_num;
        real_num = bunch.real_num;
        local_num = bunch.local_num;
	bucket_index=bunch.bucket_index;
        local_particles = new MArray2d(*(bunch.local_particles));
        state = bunch.state;
        z_period_length=bunch.z_period_length;
        z_periodic=bunch.z_periodic;
        if (bunch.converter_ptr == &(bunch.default_converter)) {
            converter_ptr = &default_converter;
        } else {
            converter_ptr = bunch.converter_ptr;
        }
    }
    return *this;
}

void
Bunch::set_particle_charge(int particle_charge)
{
    this->particle_charge = particle_charge;
}

void
Bunch::set_real_num(double real_num)
{
    this->real_num = real_num;
}

void
Bunch::set_local_num(int local_num)
{
    if (local_num > this->local_num) {
        throw std::runtime_error(
                "set_local_num can only decrease the number of particles");
    }
    this->local_num = local_num;
}

void
Bunch::update_total_num()
{
    int old_total_num = total_num;
    MPI_Allreduce(&local_num, &total_num, 1, MPI_INT, MPI_SUM,
            comm_sptr->get());
    if (old_total_num != 0.0) {
        real_num = (total_num * real_num) / old_total_num;
    } else {
        real_num = 0.0;
    }
}

void
Bunch::set_sort_period(int period)
{
    sort_period = period;
    sort_counter = period;
}

void
Bunch::sort(int index)
{
    if (index == 0) {
        Sortable2d<double*, 7, 0 > sortable(local_particles->origin(),
                local_num);
        std::sort(sortable.begin(), sortable.end(),
                Sortable2d<double*, 7, 0 >::Less());
    } else if (index == 1) {
        Sortable2d<double*, 7, 1 > sortable(local_particles->origin(),
                local_num);
        std::sort(sortable.begin(), sortable.end(),
                Sortable2d<double*, 7, 1 >::Less());
    } else if (index == 2) {
        Sortable2d<double*, 7, 2 > sortable(local_particles->origin(),
                local_num);
        std::sort(sortable.begin(), sortable.end(),
                Sortable2d<double*, 7, 2 >::Less());
    } else if (index == 3) {
        Sortable2d<double*, 7, 3 > sortable(local_particles->origin(),
                local_num);
        std::sort(sortable.begin(), sortable.end(),
                Sortable2d<double*, 7, 3 >::Less());
    } else if (index == 4) {
        Sortable2d<double*, 7, 4 > sortable(local_particles->origin(),
                local_num);
        std::sort(sortable.begin(), sortable.end(),
                Sortable2d<double*, 7, 4 >::Less());
    } else if (index == 5) {
        Sortable2d<double*, 7, 5 > sortable(local_particles->origin(),
                local_num);
        std::sort(sortable.begin(), sortable.end(),
                Sortable2d<double*, 7, 5 >::Less());
    } else {
        throw std::runtime_error("Bunch::sort: invalid index");
    }

    sort_counter = sort_period;
}

void
Bunch::periodic_sort(int index)
{
    if (sort_counter == 0) {
        sort(index);
    } else {
        --sort_counter;
    }
}

void
Bunch::set_converter(Fixed_t_z_converter &converter)
{
    this->converter_ptr = &converter;
}

void
Bunch::convert_to_state(State state)
{
    if (this->state != state) {
        if (this->state == fixed_z_lab) {
            if (state == fixed_t_lab) {
                converter_ptr->from_z_lab_to_t_lab(*this);
            }
            else if ( state == fixed_t_bunch) {
                converter_ptr->from_z_lab_to_t_bunch(*this);
            }
            // else if ( state == fixed_z_bunch) {
            //    converter_ptr->from_z_lab_to_z_bunch(*this);
            //}
            else {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 1");
            }
        }
        else if (this->state == fixed_z_bunch) {
            throw std::runtime_error("state z_bunch not implemented yet in Bunch::convert_to_state");
        }
        else if (this->state == fixed_t_lab) {
            if (state == fixed_z_lab ) {
                converter_ptr->from_t_lab_to_z_lab(*this);
            }
            //else if (state == fixed_z_bunch) {
            //    converter_ptr->from_t_lab_to_z_bunch(*this);
            //}
            else if (state == fixed_t_bunch) {
                converter_ptr->from_t_lab_to_t_bunch(*this);
            }
            else {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 2");
            }
        }
        else if (this->state == fixed_t_bunch) {
            if (state == fixed_z_lab ) {
                converter_ptr->from_t_bunch_to_z_lab(*this);
            }
            //else if (state == fixed_z_bunch ) {
            //    converter_ptr->from_t_bunch_to_z_bunch(*this);
            //}
            else if (state == fixed_t_lab ) {
                converter_ptr->from_t_bunch_to_t_lab(*this);
            }
            else {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 3");
            }
        }
    this->state =state;
    }

}



Reference_particle &
Bunch::get_reference_particle()
{
    return reference_particle;
}

Reference_particle const&
Bunch::get_reference_particle() const
{
    return reference_particle;
}

MArray2d_ref
Bunch::get_local_particles()
{
    return *local_particles;
}

Const_MArray2d_ref
Bunch::get_local_particles() const
{
    return *local_particles;
}

int
Bunch::get_particle_charge() const
{
    return particle_charge;
}

double
Bunch::get_mass() const
{
    return reference_particle.get_four_momentum().get_mass();
}



double
Bunch::get_real_num() const
{
    return real_num;
}


double
 Bunch::get_z_period_length() const
{
    return z_period_length;
}

bool
 Bunch::is_z_periodic() const
{
    return z_periodic;
}

int
Bunch::get_local_num() const
{
    return local_num;
}

int
Bunch::get_total_num() const
{
    return total_num;
}

int
Bunch::get_sort_period() const
{
    return sort_period;
}



void
Bunch::set_bucket_index(int index)
{
    this->bucket_index=index;
}

int
Bunch::get_bucket_index() const
{
return bucket_index;
}



Bunch::State
Bunch::get_state() const
{
    return state;
}

Commxx const&
Bunch::get_comm() const
{
    return *comm_sptr;
}

Commxx_sptr
Bunch::get_comm_sptr() const
{
    return comm_sptr;
}

void
Bunch::inject(Bunch const& bunch)
{
    const double weight_tolerance = 1.0e-10;
    if (std::abs(real_num / total_num - bunch.get_real_num()
            / bunch.get_total_num()) > weight_tolerance) {
        throw std::runtime_error(
                "Bunch.inject: macroparticle weight of injected bunch does not match.");
    }
    int old_local_num = local_num;
    local_num += bunch.get_local_num();
    Const_MArray2d_ref injected_particles(bunch.get_local_particles());
    MArray1d ref_state_diff(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        ref_state_diff[i] = bunch.get_reference_particle().get_state()[i]
                - reference_particle.get_state()[i];
    }
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            (*local_particles)[old_local_num + part][i]
                    = injected_particles[part][i] + ref_state_diff[i];
        }
        (*local_particles)[old_local_num + part][Bunch::id]
                = injected_particles[part][Bunch::id];
    }
    update_total_num();
}

void Bunch::check_pz2_positive()
{
    if (this->state == fixed_z_lab) {
        int local_num = get_local_num();
        MArray2d_ref particles = get_local_particles();
        for (int part = 0; part < local_num; ++part) {
            double  pzop2=(1.+particles[part][5])*(1.+particles[part][5])-
                particles[part][1]*particles[part][1]-particles[part][3]*particles[part][3];
            if (pzop2<0.)  {
                std::cout<<"pzop^2="<<pzop2<<std::endl;
                throw std::runtime_error( " check pz2:  pz square cannot be negative!");
            }

        }
    }
}


template<class Archive>
    void
    Bunch::save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP(z_period_length)
                << BOOST_SERIALIZATION_NVP(z_periodic)
                << BOOST_SERIALIZATION_NVP(reference_particle)
                << BOOST_SERIALIZATION_NVP(particle_charge)
                << BOOST_SERIALIZATION_NVP(total_num)
                << BOOST_SERIALIZATION_NVP(real_num)
                << BOOST_SERIALIZATION_NVP(bucket_index)
                << BOOST_SERIALIZATION_NVP(sort_period)
                << BOOST_SERIALIZATION_NVP(sort_counter)
                << BOOST_SERIALIZATION_NVP(state)
                << BOOST_SERIALIZATION_NVP(comm_sptr)
                << BOOST_SERIALIZATION_NVP(default_converter)
                << BOOST_SERIALIZATION_NVP(converter_ptr);
        if (comm_sptr->has_this_rank()) {
            int attempts=0;
            bool fail=true;
            while ((attempts<5) && fail){

                try {
boost::filesystem::remove(get_local_particles_serialization_path());
                    Hdf5_file file(get_local_particles_serialization_path(),
                        Hdf5_file::truncate);
                    file.write(local_num, "local_num");
                    file.write(*local_particles, "local_particles");
                    file.close();
                    fail=false;
                }
                catch(H5::Exception& he) {
                    ++attempts;
                    fail=true;
                    std::cout<<"bunch.cc: H5 Exception thrown, attempts number="
                        <<attempts<<" on rank="<<Commxx().get_rank()<<std::endl;
                    sleep(3);
                }
            }
        }
    }


template<class Archive>
    void
    Bunch::load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP(z_period_length)
                >> BOOST_SERIALIZATION_NVP(z_periodic)
                >> BOOST_SERIALIZATION_NVP(reference_particle)
                >> BOOST_SERIALIZATION_NVP(particle_charge)
                >> BOOST_SERIALIZATION_NVP(total_num)
                >> BOOST_SERIALIZATION_NVP(real_num)
                >> BOOST_SERIALIZATION_NVP(bucket_index)
                >> BOOST_SERIALIZATION_NVP(sort_period)
                >> BOOST_SERIALIZATION_NVP(sort_counter)
                >> BOOST_SERIALIZATION_NVP(state)
                >> BOOST_SERIALIZATION_NVP(comm_sptr)
                >> BOOST_SERIALIZATION_NVP(default_converter)
                >> BOOST_SERIALIZATION_NVP(converter_ptr);
        if (comm_sptr->has_this_rank()) {
            Hdf5_file file(get_local_particles_serialization_path(),
                    Hdf5_file::read_only);
            local_num = file.read<int > ("local_num");
            local_particles
                    = new MArray2d(file.read<MArray2d > ("local_particles"));
        } else {
            local_num = 0;
            local_particles = new MArray2d(boost::extents[local_num][7]);
        }
    }

template
void
Bunch::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Bunch::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Bunch::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Bunch::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Bunch::~Bunch()
{
    delete local_particles;
}
