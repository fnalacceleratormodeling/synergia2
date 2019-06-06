#include "bunch.h"
#include "synergia/utils/parallel_utils.h"
#include <boost/filesystem.hpp>
#include <boost/align/aligned_alloc.hpp>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

// particle padding based on GSVector settings
#if defined(GSV_SSE)
  const int Bunch::particle_alignment = 4;
#elif defined(GSV_AVX)
  const int Bunch::particle_alignment = 4;
#elif defined(GSV_AVX512)
  const int Bunch::particle_alignment = 8;
#elif defined(GSV_QPX)
  const int Bunch::particle_alignment = 4;
#else
  const int Bunch::particle_alignment = 4;
#endif

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
#if 0
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
#endif
}

void
Bunch::assign_spectator_ids(int local_offset)
{
#if 0

#if 0
    int global_offset, request_num;
    if (comm_sptr->get_rank() == 0) {
        request_num = total_s_num;
    } else {
        request_num = 0;
    }
    global_offset = particle_id_offset.get(request_num, *comm_sptr);
#endif

    int global_offset = 0;
    for (int i = 0; i < local_s_num; ++i) {
        (*local_s_particles)[i][id] = i + local_offset + global_offset;
    }
#endif
}



#if 0
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
#endif

std::string
Bunch::get_local_particles_serialization_path() const
{
    std::stringstream sstream;
    sstream << "local_particles_";
    sstream << bucket_index;
    sstream << ".h5";
    return get_serialization_path(sstream.str());
}

int
Bunch::calculate_aligned_pos(int num)
{
    if (particle_alignment <= 0)
    {
        throw std::runtime_error("Bunch::calculate_aligned_pos() invalid particle_alignment value");
    }

    if (num < 0)
    {
        throw std::runtime_error("Bunch::calculate_aligned_pos() invalid num value");
    }

    if (num == 0)
    {
        return 0;
    }
    else
    {
        if (num % particle_alignment == 0) 
        {
            return num;
        } 
        else 
        {
            return num + particle_alignment - (num % particle_alignment);
        }
    }
}

int
Bunch::calculate_padding_size(int num)
{
    if (particle_alignment <= 0)
    {
        throw std::runtime_error("Bunch::calculate_padding_size() invalid particle_alignment value");
    }

    if (num < 0)
    {
        throw std::runtime_error("Bunch::calculate_padding_size() invalid num value");
    }

    if (num == 0)
    {
        return 0;
    }
    else
    {
        if (num % particle_alignment == 0) 
        {
            return particle_alignment;
        } 
        else 
        {
            return particle_alignment * 2 - num % particle_alignment;
        }
    }
}

void
Bunch::construct(int total_num, double real_num, int total_s_num)
{
    sort_counter = 0;
    sort_period = 10000;

    this->total_num = total_num;
    this->real_num = real_num;
    this->total_s_num = total_s_num;

    //state = fixed_z_lab;
    //converter_ptr = &default_converter;

    if (comm_sptr->has_this_rank()) 
    {
        // real particles
        // ----------------------------------------------------------------------
        std::vector<int> offsets(comm_sptr->get_size());
        std::vector<int> counts(comm_sptr->get_size());
        decompose_1d(*comm_sptr, total_num, offsets, counts);

        local_num         = counts[comm_sptr->get_rank()];
        local_num_aligned = calculate_aligned_pos(local_num);
        local_num_padded  = local_num + calculate_padding_size(local_num);
        local_num_slots   = local_num_padded;

        // allocate
        Kokkos::resize(parts, local_num_slots);
        hparts = Kokkos::create_mirror_view(parts);

#if 0
        // zero
        #pragma omp parallel for
        for (int i=0; i<local_num_slots; ++i)
        {
            for(int j=0; j<7; ++j)
            {
                (*local_particles)[i][j] = 0.0;
            }
        }
#endif

        // id
        assign_ids(offsets[comm_sptr->get_rank()]);

        // spectator particles
        // ----------------------------------------------------------------------
        std::vector<int> s_offsets(comm_sptr->get_size());
        std::vector<int> s_counts(comm_sptr->get_size());
        decompose_1d(*comm_sptr, total_s_num, s_offsets, s_counts);

        local_s_num         = s_counts[comm_sptr->get_rank()];
        local_s_num_aligned = calculate_aligned_pos(local_s_num);
        local_s_num_padded  = local_s_num + calculate_padding_size(local_s_num);
        local_s_num_slots   = local_s_num_padded;

        // allocate
        Kokkos::resize(sparts, local_s_num_slots);
        hsparts = Kokkos::create_mirror_view(sparts);

#if 0
        // reset
        #pragma omp parallel for
        for (int i=0; i<local_s_num_slots; ++i)
        {
            for(int j=0; j<7; ++j)
            {
                (*local_s_particles)[i][j] = 0.0;
            }
        }
#endif

        // id
        assign_spectator_ids(s_offsets[comm_sptr->get_rank()]);
    } 
    else 
    {
        local_num = 0;
        local_num_aligned = 0;
        local_num_padded = 0;
        local_num_slots = 0;

        local_s_num = 0;
        local_s_num_aligned = 0;
        local_s_num_padded = 0;
        local_s_num_slots = 0;
    }
}

Bunch::Bunch(
        Reference_particle const& reference_particle, 
        int total_num, 
        double real_num, 
        Commxx_sptr comm_sptr) 
    : longitudinal_extent(0.0)
    , z_periodic(false)
    , longitudinal_aperture(false)
    , reference_particle(reference_particle)
    , design_reference_particle(reference_particle)
    , particle_charge(reference_particle.get_charge())

    , local_num(0)
    , local_num_aligned(0)
    , local_num_padded(0)
    , local_num_slots(0)

    , local_s_num(0)
    , local_s_num_aligned(0)
    , local_s_num_padded(0)
    , local_s_num_slots(0)

    , total_num(total_num)
    , total_s_num(0)

    , real_num(real_num)

    , parts("particles", local_num_slots)
    , hparts(Kokkos::create_mirror_view(parts))

    , sparts("spectator particles", local_s_num_slots)
    , hsparts(Kokkos::create_mirror_view(sparts))

    , bucket_index(0)
    , bucket_index_assigned(false)
    , sort_period(10000)
    , sort_counter(0)
    , comm_sptr(comm_sptr)
{
    construct(total_num, real_num, 0);
}

Bunch::Bunch(
        Reference_particle const& reference_particle, 
        int total_num, 
        int total_spectator_num, 
        double real_num, 
        Commxx_sptr comm_sptr) 
    : longitudinal_extent(0.0)
    , z_periodic(false)
    , longitudinal_aperture(false)
    , reference_particle(reference_particle)
    , design_reference_particle(reference_particle)
    , particle_charge(reference_particle.get_charge())

    , local_num(0)
    , local_num_aligned(0)
    , local_num_padded(0)
    , local_num_slots(0)

    , local_s_num(0)
    , local_s_num_aligned(0)
    , local_s_num_padded(0)
    , local_s_num_slots(0)

    , total_num(total_num)
    , total_s_num(total_spectator_num)

    , real_num(real_num)

    , parts("particles", local_num_slots)
    , hparts(Kokkos::create_mirror_view(parts))

    , sparts("spectator particles", local_s_num_slots)
    , hsparts(Kokkos::create_mirror_view(sparts))

    , bucket_index(0)
    , bucket_index_assigned(false)
    , sort_period(10000)
    , sort_counter(0)
    , comm_sptr(comm_sptr)
{
    construct(total_num, real_num, total_spectator_num);
}

Bunch::Bunch() 
    : longitudinal_extent(0.0)
    , z_periodic(false)
    , longitudinal_aperture(false)
    , reference_particle()
    , design_reference_particle()
    , particle_charge(1)

    , local_num(0)
    , local_num_aligned(0)
    , local_num_padded(0)
    , local_num_slots(0)

    , local_s_num(0)
    , local_s_num_aligned(0)
    , local_s_num_padded(0)
    , local_s_num_slots(0)

    , total_num(0)
    , total_s_num(0)

    , real_num(0)

    , parts("particles", local_num_slots)
    , hparts(Kokkos::create_mirror_view(parts))

    , sparts("spectator particles", local_s_num_slots)
    , hsparts(Kokkos::create_mirror_view(sparts))

    , bucket_index(0)
    , bucket_index_assigned(false)
    , sort_period(10000)
    , sort_counter(0)
    , comm_sptr()
{
}

#if 0
Bunch &
Bunch::operator=(Bunch const& bunch)
{
    if (this != &bunch) 
    {
        reference_particle = bunch.reference_particle;
        design_reference_particle = bunch.design_reference_particle;

        comm_sptr = bunch.comm_sptr;

        particle_charge = bunch.particle_charge;

        local_num = bunch.local_num;
        local_num_aligned = bunch.local_num_aligned;
        local_num_padded = bunch.local_num_padded;
        local_num_slots = bunch.local_num_slots;

        local_s_num = bunch.local_s_num;
        local_s_num_aligned = bunch.local_s_num_aligned;
        local_s_num_padded = bunch.local_s_num_padded;
        local_s_num_slots = bunch.local_s_num_slots;

        total_num = bunch.total_num;
        total_s_num = bunch.total_s_num;

        real_num = bunch.real_num;

        bucket_index = bunch.bucket_index;
        bucket_index_assigned = bunch.bucket_index_assigned;

        // delete current storages
        if (storage) boost::alignment::aligned_free(storage);
        if (local_particles) delete local_particles;

        if (s_storage) boost::alignment::aligned_free(s_storage);
        if (local_s_particles) delete local_s_particles;

        // and allocate new
        storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_num_slots * 7 * sizeof(double));
        memcpy(storage, bunch.storage, sizeof(double)*local_num_slots*7);
        local_particles = new MArray2d_ref(storage, boost::extents[local_num_slots][7], boost::fortran_storage_order());

        s_storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_s_num_slots * 7 * sizeof(double));
        memcpy(s_storage, bunch.s_storage, sizeof(double)*local_s_num_slots *7);
        local_s_particles = new MArray2d_ref(storage, boost::extents[local_s_num_slots][7], boost::fortran_storage_order());

        longitudinal_extent = bunch.longitudinal_extent;
        z_periodic = bunch.z_periodic;
        longitudinal_aperture = bunch.longitudinal_aperture;
    }

    return *this;
}
#endif

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
Bunch::expand_local_num(int num, int added_lost)
{
#if 0
    // keep the previous values
    int prev_local_num = local_num;
    int prev_local_num_padded = local_num_padded;

    int local_num_lost = local_num_slots - local_num_padded;
    int total_num_lost = local_num_lost + added_lost;

    double * prev_storage = storage;
    MArray2d_ref * prev_local_particles = local_particles;

    // update the pointers
    local_num = num;
    local_num_aligned = calculate_aligned_pos(local_num);
    local_num_padded = local_num + calculate_padding_size(local_num + total_num_lost);
    local_num_slots = local_num_padded + total_num_lost;

    // allocate for new storage
    storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_num_slots * 7 * sizeof(double));
    local_particles = new MArray2d_ref(storage, boost::extents[local_num_slots][7], boost::fortran_storage_order());

    // copy the particle data over
    for (int i = 0; i < prev_local_num; ++i) 
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[i][j] = (*prev_local_particles)[i][j];
        }
    }

    // set the coordinates of extended and padding particles to 0
    // TODO: what should be the id for the extended particles
    for (int i = prev_local_num; i < local_num_padded; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[i][j] = 0.0;
        }
    }

    // copy over lost particles
    for (int i=0; i<local_num_lost; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[local_num_padded + i][j] = 
                (*prev_local_particles)[prev_local_num_padded + i][j];
        }
    }

    // set additional lost particle data to 0
    for (int i = local_num_padded + local_num_lost; i < local_num_slots; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[i][j] = 0.0;
        }
    }

    if (prev_storage) boost::alignment::aligned_free(prev_storage);
    if (prev_local_particles) delete prev_local_particles;
#endif
}

void
Bunch::expand_local_spectator_num(int s_num, int added_lost)
{
#if 0
    // keep the previous values
    int prev_local_s_num = local_s_num;
    int prev_local_s_num_padded = local_s_num_padded;

    int local_s_num_lost = local_s_num_slots - local_s_num_padded;
    int total_s_num_lost = local_s_num_lost + added_lost;

    double * prev_s_storage = s_storage;
    MArray2d_ref * prev_local_s_particles = local_s_particles;

    // update the pointers
    local_s_num = s_num;
    local_s_num_aligned = calculate_aligned_pos(local_s_num);
    local_s_num_padded = local_s_num + calculate_padding_size(local_s_num + total_s_num_lost);
    local_s_num_slots = local_s_num_padded + total_s_num_lost;

    // allocate for new storage
    s_storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_s_num_slots * 7 * sizeof(double));
    local_s_particles = new MArray2d_ref(s_storage, boost::extents[local_s_num_slots][7], boost::fortran_storage_order());

    // copy the particle data over
    for (int i = 0; i < prev_local_s_num; ++i) 
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_s_particles)[i][j] = (*prev_local_s_particles)[i][j];
        }
    }

    // set the coordinates of extended and padding particles to 0
    // TODO: what should be the id for the extended particles
    for (int i = prev_local_s_num; i < local_s_num_padded; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_s_particles)[i][j] = 0.0;
        }
    }

    // copy over lost particles
    for (int i=0; i<local_s_num_lost; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_s_particles)[local_s_num_padded + i][j] = 
                (*prev_local_s_particles)[prev_local_s_num_padded + i][j];
        }
    }

    // set additional lost particle data to 0
    for (int i = local_s_num_padded + local_s_num_lost; i < local_s_num_slots; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_s_particles)[i][j] = 0.0;
        }
    }

    if (prev_s_storage) boost::alignment::aligned_free(prev_s_storage);
    if (prev_local_s_particles) delete prev_local_s_particles;
#endif
}

void
Bunch::set_local_num(int num)
{
#if 0
    // make sure the new local_num is never less than 0
    if (num < 0) num = 0;

    // re-allocate depending on the new size
    if (num <= local_num_aligned)
    {
        // previous values
        int prev_local_num = local_num;

        // no need to resize the array, only move the pointers
        local_num = num;

        // update local_num_aligned
        local_num_aligned = calculate_aligned_pos(local_num);

        // clear the particle data (from local_num to prev_local_num)
        // note this only happens when the new local_num is smaller than the old one
        for (int i = local_num; i < prev_local_num; ++i)
        {
            for (int j=0; j<7; ++j)
            {
                (*local_particles)[i][j] = 0.0;
            }
        }
    }
    else
    {
        // expand the local particle array, no additional lost particle slots
        expand_local_num(num, 0);
    }
#endif
}

void
Bunch::set_local_spectator_num(int s_num)
{
#if 0
    // make sure the new local_s_num is never less than 0
    if (s_num < 0) s_num = 0;

    // re-allocate depending on the new size
    if (s_num <= local_s_num_aligned)
    {
        // previous values
        int prev_local_s_num = local_s_num;

        // no need to resize the array, only move the pointers
        local_s_num = s_num;

        // update local_s_num_aligned
        local_s_num_aligned = calculate_aligned_pos(local_s_num);

        // clear the particle data (from local_s_num to prev_local_s_num)
        // note this only happens when the new local_s_num is smaller than the old one
        for (int i = local_s_num; i < prev_local_s_num; ++i)
        {
            for (int j=0; j<7; ++j)
            {
                (*local_s_particles)[i][j] = 0.0;
            }
        }
    }
    else
    {
        // expand the local spectator particle array, no additional lost particle slots
        expand_local_spectator_num(s_num, 0);
    }
#endif
}

void
Bunch::update_total_num()
{
#if 0
    // total real particle number
    int old_total_num = total_num;
    MPI_Allreduce(&local_num, &total_num, 1, MPI_INT, MPI_SUM, comm_sptr->get());

    // total spectator particle number
    MPI_Allreduce(&local_s_num, &total_s_num, 1, MPI_INT, MPI_SUM, comm_sptr->get());

    if (old_total_num != 0) 
    {
        real_num = (total_num * real_num) / old_total_num;
    } 
    else 
    {
        real_num = 0.0;
    }
#endif
}

void
Bunch::set_total_num(int totalnum)
{
#if 0
    int old_total_num = total_num;
    total_num = totalnum;

    if (old_total_num != 0) 
    {
        real_num = (total_num * real_num) / old_total_num;
    } 
    else 
    {
        real_num = 0.0;
    }
#endif
}

#if 0
namespace {
double * semi_global_t;
size_t semi_global_start_pos;

inline bool do_compare(unsigned int const& a, unsigned int const& b)
{
    bool retval = semi_global_t[semi_global_start_pos+a] <
            semi_global_t[semi_global_start_pos+b];
    return retval;
}

void do_sort(double * t, size_t rows, size_t cols, size_t cols_padded, size_t ord_col)
{
    semi_global_t = t;
    std::vector<unsigned int> index(cols);
    // c++ 11
    // unsigned int ind=0;
    //generate(index.begin(),index.end(), [&]() { return ind++; });
    for(int i=0; i<cols; ++i) {
        index[i] = i;
    }
    semi_global_start_pos = ord_col*cols_padded;
    std::sort(index.begin(),index.end(), &do_compare);

    // swap all values in each row according to index order
    for(size_t r=0; r<rows; ++r) {
        double *start = t+(r*cols_padded), *end = t+((r+1)*cols_padded);
        std::vector<double> temp(start, end);
        for(size_t i=0; i<index.size(); ++i) {
            start[i] = temp[index[i]];
        }
    }
}
}   // ends namespace

void
Bunch::set_sort_period(int period)
{
    sort_period = period;
    sort_counter = period;
}

void
Bunch::sort(int index)
{
    if ((index<0) || (index>6)) {
        throw std::runtime_error("Bunch::sort: invalid index");
    }
    do_sort(local_particles->origin(), 7, local_num, local_num_slots, index);
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
#endif

#if 0
void
Bunch::set_converter(Fixed_t_z_converter &converter)
{
    //this->converter_ptr = &converter;
}
#endif

#if 0
void
Bunch::convert_to_state(State state)
{
    if (this->state != state) 
    {
        if (this->state == fixed_z_lab) 
        {
            if (state == fixed_t_lab) 
            {
                converter_ptr->from_z_lab_to_t_lab(*this);
            }
            else if ( state == fixed_t_bunch) 
            {
                converter_ptr->from_z_lab_to_t_bunch(*this);
            }
            // else if ( state == fixed_z_bunch) {
            //    converter_ptr->from_z_lab_to_z_bunch(*this);
            //}
            else 
            {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 1");
            }
        }
        else if (this->state == fixed_z_bunch) 
        {
            throw std::runtime_error("state z_bunch not implemented yet in Bunch::convert_to_state");
        }
        else if (this->state == fixed_t_lab) 
        {
            if (state == fixed_z_lab ) 
            {
                converter_ptr->from_t_lab_to_z_lab(*this);
            }
            //else if (state == fixed_z_bunch) {
            //    converter_ptr->from_t_lab_to_z_bunch(*this);
            //}
            else if (state == fixed_t_bunch) 
            {
                converter_ptr->from_t_lab_to_t_bunch(*this);
            }
            else 
            {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 2");
            }
        }
        else if (this->state == fixed_t_bunch) 
        {
            if (state == fixed_z_lab ) 
            {
                converter_ptr->from_t_bunch_to_z_lab(*this);
            }
            //else if (state == fixed_z_bunch ) {
            //    converter_ptr->from_t_bunch_to_z_bunch(*this);
            //}
            else if (state == fixed_t_lab ) 
            {
                converter_ptr->from_t_bunch_to_t_lab(*this);
            }
            else 
            {
                std::cout<<" state to convert to="<<state<<std::endl;
                std::cout<<" initial state ="<<this->state<<std::endl;
                throw std::runtime_error("Unknown state in Bunch::convert_to_state, case 3");
            }
        }

        this->state = state;
    }
}
#endif



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

Reference_particle &
Bunch::get_design_reference_particle()
{
    return design_reference_particle;
}

Reference_particle const&
Bunch::get_design_reference_particle() const
{
    return design_reference_particle;
}

void
Bunch::set_design_reference_particle(Reference_particle const & ref_part)
{
    design_reference_particle = ref_part;
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
    return longitudinal_extent;
}

void
 Bunch::set_z_period_length(double z_period_length)
{
    if (longitudinal_aperture)  throw std::runtime_error("longitudinal_aperture is true, cannot make the bunch z periodic");
    this->longitudinal_extent=z_period_length;
    this->z_periodic=true;
}

bool
 Bunch::is_z_periodic() const
{
    return z_periodic;
}

double
Bunch::get_longitudinal_aperture_length() const
{
    return longitudinal_extent;
}

void
Bunch::set_longitudinal_aperture_length(double longitudinal_extent)
{
    if (z_periodic)  throw std::runtime_error("z_periodic is true, cannot put a longitudinal_aperture");
    this->longitudinal_extent=longitudinal_extent;
    this->longitudinal_aperture=true;
}

bool
Bunch::has_longitudinal_aperture() const
{
  return longitudinal_aperture;
}

int
Bunch::get_local_num() const
{
    return local_num;
}

int
Bunch::get_local_num_aligned() const
{
    return local_num_aligned;
}

int
Bunch::get_local_num_padded() const
{
    return local_num_padded;
}

int
Bunch::get_local_num_padding() const
{
    return local_num_padded - local_num;
}

int
Bunch::get_local_num_lost() const
{
    return local_num_slots - local_num_padded;
}

int
Bunch::get_local_num_slots() const
{
    return local_num_slots;
}

int
Bunch::get_total_num() const
{
    return total_num;
}

int
Bunch::get_local_spectator_num() const
{
    return local_s_num;
}

int
Bunch::get_local_spectator_num_aligned() const
{
    return local_s_num_aligned;
}

int
Bunch::get_local_spectator_num_padded() const
{
    return local_s_num_padded;
}

int
Bunch::get_local_spectator_num_padding() const
{
    return local_s_num_padded - local_s_num;
}

int
Bunch::get_local_spectator_num_lost() const
{
    return local_s_num_slots - local_s_num_padded;
}

int
Bunch::get_local_spectator_num_slots() const
{
    return local_s_num_slots;
}

int
Bunch::get_total_spectator_num() const
{
    return total_s_num;
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
    this->bucket_index_assigned=true;
}

int
Bunch::get_bucket_index() const
{
  if (!bucket_index_assigned)  throw std::runtime_error("bucket index has not been assigned yet");
  return bucket_index;
}

bool
Bunch::is_bucket_index_assigned() const
{
  return bucket_index_assigned;
}


Bunch::State
Bunch::get_state() const
{
    return fixed_z_lab;
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
    throw std::runtime_error("Bunch::inject not implemented");

#if 0
    const double weight_tolerance = 1.0e-14;
    const double particle_tolerance = 1.0e-14;

    // The charge and mass of the bunch particles must match
    if (particle_charge != bunch.get_particle_charge()) 
    {
        throw std::runtime_error(
                "Bunch.inject: bunch particle charges do not match.");
    }

    if (std::abs(reference_particle.get_four_momentum().get_mass()/
                 bunch.get_reference_particle().get_four_momentum().get_mass() - 1.0) > particle_tolerance) 
    {
        throw std::runtime_error(
                "Bunch:inject: bunch particle masses do not match.");
    }

    // can only check particle weight if total_num is nonzero
    if (total_num == 0) 
    {
        // target bunch is empty.  Set the weights from the injected bunch
        real_num = bunch.get_real_num();
        total_num = bunch.get_total_num();
    } 
    else 
    {
        double wgt1 = real_num/total_num;
        double wgt2 = bunch.get_real_num()/bunch.get_total_num();

        if (std::abs(wgt1/wgt2 - 1.0) > weight_tolerance) 
        {
            throw std::runtime_error(
                "Bunch.inject: macroparticle weight of injected bunch does not match.");
        }
    }

    ConstParticles injected_particles(bunch.get_local_particles());
    ConstParticles injected_spectator_particles(bunch.get_local_spectator_particles());

    double target_momentum = reference_particle.get_momentum();
    double injected_momentum = bunch.get_reference_particle().get_momentum();

    MArray1d ref_state_diff(boost::extents[6]);
    MArray1d target_state(boost::extents[6]);
    MArray1d injected_state(boost::extents[6]);

    for (int i = 0; i < 6; ++i) 
    {
        ref_state_diff[i] = bunch.get_reference_particle().get_state()[i]
                - reference_particle.get_state()[i];
    }

    for (int i = 0; i < 6; ++i) 
    {
        target_state[i] = reference_particle.get_state()[i];
        injected_state[i] = bunch.get_reference_particle().get_state()[i];
    }

    // real particles
    int old_local_num = local_num;
    int old_local_num_lost = get_local_num_lost();

    // reallocate the particle array
    expand_local_num(
            old_local_num + bunch.get_local_num(), 
            bunch.get_local_num_lost() );

    for (int part = 0; part < bunch.get_local_num(); ++part) 
    {
        // space-like coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_particles)[old_local_num + part][i]
                    = injected_particles[part][i] + ref_state_diff[i];
        }

        // npx and npy coordinates are scaled with p_ref which can be different
        // for different bunches
        for (int i = 1; i < 4; i += 2) 
        {
            (*local_particles)[old_local_num + part][i] =
                    (injected_momentum/target_momentum) *
                    (injected_particles[part][i] - injected_state[i]) + target_state[i];
        }

        // ndp coordinate is delta-p scaled with pref
        (*local_particles)[old_local_num + part][5] =
                (injected_momentum/target_momentum) *
                (1.0 + injected_particles[part][5] - injected_state[5]) + target_state[5] - 1.0;

        (*local_particles)[old_local_num + part][Bunch::id]
                = injected_particles[part][Bunch::id];
    }

    // copy the lost particles from bunch
    for (int p = 0; p < bunch.get_local_num_lost(); ++p)
    {
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_particles)[local_num_padded + old_local_num_lost + p][i] =
                injected_particles[bunch.get_local_num_padded() + p][i];
        }
    }

    // spectator particles
    int old_local_s_num = local_s_num;
    int old_local_s_num_lost = get_local_spectator_num_lost();

    // reallocate the spectator particle array
    expand_local_spectator_num(
            old_local_s_num + bunch.get_local_spectator_num(), 
            bunch.get_local_spectator_num_lost() );

    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        // space-like coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_s_particles)[old_local_s_num + part][i]
                    = injected_spectator_particles[part][i] + ref_state_diff[i];
        }

        // npx and npy coordinates are scaled with p_ref which can be different
        // for different bunches
        for (int i = 1; i < 4; i += 2) 
        {
            (*local_s_particles)[old_local_s_num + part][i] =
                    (injected_momentum/target_momentum) *
                    (injected_spectator_particles[part][i] - injected_state[i]) + target_state[i];
        }

        // ndp coordinate is delta-p scaled with pref
        (*local_s_particles)[old_local_s_num + part][5] =
                (injected_momentum/target_momentum) *
                (1.0 + injected_spectator_particles[part][5] - injected_state[5]) + target_state[5] - 1.0;

        (*local_s_particles)[old_local_s_num + part][Bunch::id]
                = injected_spectator_particles[part][Bunch::id];
    }

    // copy the lost spectator particles from bunch
    for (int p = 0; p < bunch.get_local_spectator_num_lost(); ++p)
    {
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_s_particles)[local_s_num_padded + old_local_s_num_lost + p][i] =
                injected_spectator_particles[bunch.get_local_spectator_num_padded() + p][i];
        }
    }

    // update total number, for both real and spectator particles
    update_total_num();
#endif
}

void
Bunch::read_file(std::string const & filename)
{
    throw std::runtime_error("Bunch::read_file() no implemented");

#if 0
   if (comm_sptr->has_this_rank()) 
   {
        Hdf5_file file(filename, Hdf5_file::read_only);
        MArray2d* read_particles = new MArray2d(file.read<MArray2d>("particles"));

        int num_particles = read_particles->shape()[0];
        if (total_num != num_particles) 
        {
            throw std::runtime_error( 
                    " the initial bunch file has a different number of particles");
        }

        std::vector<int> offsets(comm_sptr->get_size());
        std::vector<int> counts(comm_sptr->get_size());
        decompose_1d(*comm_sptr, total_num, offsets, counts);

        if (local_num !=  counts[comm_sptr->get_rank()]) 
        {
            throw std::runtime_error( 
                    " local_num incompatibility when initializing the bunch");
        }

        int offset = offsets[comm_sptr->get_rank()];

        for (int part = 0; part < local_num; ++part) 
        {
            int rpart = part + offset;

            for (int i = 0; i < 7; ++i) 
            {
                (*local_particles)[part][i]=(*read_particles)[rpart][i];
            }
        }
   }
#endif
}

void Bunch::check_pz2_positive()
{
    throw std::runtime_error("Bunch::check_pz2_positive() not implemented");

#if 0
    if (this->state == fixed_z_lab) 
    {
        int local_num = get_local_num();
        MArray2d_ref particles = get_local_particles();

        for (int part = 0; part < local_num; ++part) 
        {
            double  pzop2=(1.+particles[part][5])*(1.+particles[part][5])-
                particles[part][1]*particles[part][1]-particles[part][3]*particles[part][3];

            if (pzop2<0.)  
            {
                std::cout<<"pzop^2="<<pzop2<<std::endl;
                throw std::runtime_error( " check pz2:  pz square cannot be negative!");
            }
        }
    }
#endif
}

#if 0
template<class Archive>
void
Bunch::save(Archive & ar, const unsigned int version) const
{
    ar << CEREAL_NVP(longitudinal_extent)
       << CEREAL_NVP(z_periodic)
       << CEREAL_NVP(longitudinal_aperture)

       << CEREAL_NVP(reference_particle)
       << CEREAL_NVP(design_reference_particle)
       << CEREAL_NVP(particle_charge)

       << CEREAL_NVP(total_num)
       << CEREAL_NVP(total_s_num)
       << CEREAL_NVP(real_num)

       << CEREAL_NVP(bucket_index)
       << CEREAL_NVP(bucket_index_assigned)

       << CEREAL_NVP(sort_period)
       << CEREAL_NVP(sort_counter)

       << CEREAL_NVP(comm_sptr)

       //<< CEREAL_NVP(state)
       //<< CEREAL_NVP(default_converter)
       //<< CEREAL_NVP(converter_ptr)
       ;

    if (comm_sptr->has_this_rank()) 
    {
        int attempts = 0;
        bool fail = true;

        while ((attempts<5) && fail)
        {
            try 
            {
                boost::filesystem::remove(get_local_particles_serialization_path());
                Hdf5_file file(get_local_particles_serialization_path(), Hdf5_file::truncate);

                file.write(local_num, "local_num");
                file.write(local_num_aligned, "local_num_aligned");
                file.write(local_num_padded, "local_num_padded");
                file.write(local_num_slots, "local_num_slots");
                file.write(storage, local_num_slots*7, "local_storage");

                file.write(local_s_num, "local_s_num");
                file.write(local_s_num_aligned, "local_s_num_aligned");
                file.write(local_s_num_padded, "local_s_num_padded");
                file.write(local_s_num_slots, "local_s_num_slots");
                file.write(s_storage, local_s_num_slots*7, "local_s_storage");

                file.close();
                fail=false;
            }
            catch(Hdf5_exception & he) 
            {
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
    ar >> CEREAL_NVP(longitudinal_extent)
       >> CEREAL_NVP(z_periodic)
       >> CEREAL_NVP(longitudinal_aperture)

       >> CEREAL_NVP(reference_particle)
       >> CEREAL_NVP(design_reference_particle)
       >> CEREAL_NVP(particle_charge)

       >> CEREAL_NVP(total_num)
       >> CEREAL_NVP(total_s_num)
       >> CEREAL_NVP(real_num)

       >> CEREAL_NVP(bucket_index)
       >> CEREAL_NVP(bucket_index_assigned)

       >> CEREAL_NVP(sort_period)
       >> CEREAL_NVP(sort_counter)

       >> CEREAL_NVP(comm_sptr)

       //>> CEREAL_NVP(state)
       //>> CEREAL_NVP(default_converter)
       //>> CEREAL_NVP(converter_ptr)
       ;

    if (comm_sptr->has_this_rank()) 
    {
        Hdf5_file file(get_local_particles_serialization_path(), Hdf5_file::read_only);

        local_num = file.read<int> ("local_num");
        local_num_aligned = file.read<int> ("local_num_aligned");
        local_num_padded = file.read<int> ("local_num_padded");
        local_num_slots = file.read<int> ("local_num_slots");

        local_s_num = file.read<int> ("local_s_num");
        local_s_num_aligned = file.read<int> ("local_s_num_aligned");
        local_s_num_padded = file.read<int> ("local_s_num_padded");
        local_s_num_slots = file.read<int> ("local_s_num_slots");

        storage = file.read<double *>("local_storage");
        local_particles = new MArray2d_ref(storage, boost::extents[local_num_slots][7], boost::fortran_storage_order());

        s_storage = file.read<double *>("local_s_storage");
        local_s_particles = new MArray2d_ref(s_storage, boost::extents[local_s_num_slots][7], boost::fortran_storage_order());
    } 
    else 
    {
        local_num = 0;
        local_num_aligned = 0;
        local_num_padded = 0;
        local_num_slots = 0;

        local_s_num = 0;
        local_s_num_aligned = 0;
        local_s_num_padded = 0;
        local_s_num_slots = 0;

        storage = NULL;
        local_particles = new MArray2d_ref(storage, boost::extents[0][7], boost::fortran_storage_order());

        s_storage = NULL;
        local_s_particles = new MArray2d_ref(s_storage, boost::extents[0][7], boost::fortran_storage_order());
    }
}

template
void
Bunch::save<cereal::BinaryOutputArchive >(
        cereal::BinaryOutputArchive & ar, const unsigned int version) const;

template
void
Bunch::save<cereal::XMLOutputArchive >(
        cereal::XMLOutputArchive & ar, const unsigned int version) const;

template
void
Bunch::load<cereal::BinaryInputArchive >(
        cereal::BinaryInputArchive & ar, const unsigned int version);

template
void
Bunch::load<cereal::XMLInputArchive >(
        cereal::XMLInputArchive & ar, const unsigned int version);

#endif
