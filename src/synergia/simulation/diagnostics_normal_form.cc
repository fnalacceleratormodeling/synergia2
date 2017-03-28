#include "diagnostics_normal_form.h"

const char Diagnostics_normal_form::name[] = "diagnostics_normal_form";

Diagnostics_normal_form::Diagnostics_normal_form(Fast_normal_form_sptr fast_nf_sptr, std::string const& filename,
        std::string const& local_dir) :
        Diagnostics_normal_form::Diagnostics(Diagnostics_normal_form::name, filename,
                local_dir), fnf_sptr(fast_nf_sptr),  have_writers(false), writer_s_n(0), writer_repetition(
                0), writer_s(0), writer_num_particles(0),  aa2(boost::extents[3]), writer_aa2(0)
               
{
}

Diagnostics_normal_form::Diagnostics_normal_form(): have_writers(false)
{
}


bool
Diagnostics_normal_form::is_serial() const
{
    return true;
}



Diagnostics_normal_form::~Diagnostics_normal_form()
{
     if (have_writers) {
         delete writer_aa2;
         delete writer_num_particles;
         delete writer_s;
         delete writer_repetition;
         delete writer_s_n;
     }
}


void
Diagnostics_normal_form::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
       
        writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
        writer_repetition = new Hdf5_serial_writer<int > (file_sptr,
                "repetition");
        writer_s = new Hdf5_serial_writer<double > (file_sptr,
                "s");
        writer_num_particles = new Hdf5_serial_writer<int > (file_sptr,
                "num_particles");       
        writer_aa2 = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "aa2");
        
        have_writers = true;
    }
}

void
Diagnostics_normal_form::update_aa2()
{   
     int num_part=get_bunch().get_local_num();
     MArray2d particles_copy(boost::extents[num_part][7]);
     particles_copy=get_bunch().get_local_particles();
     fnf_sptr->convert_xyz_to_normal(particles_copy);
     
     MArray1d sum2(boost::extents[3]);
     for (int j=0;j<3;++j){ 
         sum2[j]=0.;
     }
     
    
     for (int p=0;p<num_part;++p){
         for (int j=0;j<3;++j){ 
               sum2[j] += particles_copy[p][2*j]*particles_copy[p][2*j]+particles_copy[p][2*j+1]*particles_copy[p][2*j+1];             
         }
     }
     MPI_Allreduce(sum2.origin(), aa2.origin(), 3, MPI_DOUBLE, MPI_SUM,
            get_bunch().get_comm().get());

     for (int j=0;j<3;++j){
        aa2[j]  /= get_bunch().get_total_num();
     }
     
     
}

void
Diagnostics_normal_form::update()
{   if (get_bunch().get_comm().has_this_rank()){
      get_bunch().convert_to_state(Bunch::fixed_z_lab);
      s_n = get_bunch().get_reference_particle().get_s_n();
      repetition = get_bunch().get_reference_particle().get_repetition();
      s = get_bunch().get_reference_particle().get_s();
      num_particles = get_bunch().get_total_num();
      update_aa2();      
    }
}

void
Diagnostics_normal_form::write()
{
    if (get_bunch().get_comm().has_this_rank()){
      if (get_write_helper().write_locally()) {
      init_writers(get_write_helper().get_hdf5_file_sptr());
      writer_s_n->append(s_n);
      writer_repetition->append(repetition);
      writer_s->append(s);   
      writer_num_particles->append(num_particles);
      writer_aa2->append(aa2);
      get_write_helper().finish_write();
      }
    }
}

template<class Archive>
    void
    Diagnostics_normal_form::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                & BOOST_SERIALIZATION_NVP(fnf_sptr)
                & BOOST_SERIALIZATION_NVP(have_writers)
                & BOOST_SERIALIZATION_NVP(s_n)
                & BOOST_SERIALIZATION_NVP(writer_s_n)
                & BOOST_SERIALIZATION_NVP(repetition)
                & BOOST_SERIALIZATION_NVP(writer_repetition)
                & BOOST_SERIALIZATION_NVP(s)
                & BOOST_SERIALIZATION_NVP(writer_s) 
                & BOOST_SERIALIZATION_NVP(num_particles)
                & BOOST_SERIALIZATION_NVP(writer_num_particles)
                & BOOST_SERIALIZATION_NVP(aa2)
                & BOOST_SERIALIZATION_NVP(writer_aa2);
    }

template
void
Diagnostics_normal_form::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_normal_form::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_normal_form::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_normal_form::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_normal_form)
