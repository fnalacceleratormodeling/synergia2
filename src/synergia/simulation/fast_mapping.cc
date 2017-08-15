#include "fast_mapping.h"
//#include "synergia/lattice/chef_utils.h"

//#define MANUAL_LOOP_UNROLL shockingly_yes

std::string
read_line_ignoring_comments(std::ifstream &file)
{
    std::string line;
    getline(file, line);
    if (line.size() == 0) {
        return line;
    }
    while (((line[0] == '#') || (line.size() == 0))) {
        getline(file, line);
    }
    return line;
}

Fast_mapping_term::Fast_mapping_term(int order) : i(order + 1)
{
    this->the_order = order;
}

Fast_mapping_term::Fast_mapping_term(std::ifstream & stream)
{
    std::stringstream sstream(read_line_ignoring_comments(stream));
    sstream >> the_order;
    i.resize(the_order + 1);
    for (int j = 0; j < the_order + 1; ++j) {
        sstream >> i[j];
    }
    sstream >> the_coeff;
}

Fast_mapping_term::Fast_mapping_term()
{
}

template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template<class Archive>
    void
    Fast_mapping_term::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(the_coeff);
        ar & BOOST_SERIALIZATION_NVP(the_order);
        ar & BOOST_SERIALIZATION_NVP(i);
    }

template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Fast_mapping_term::~Fast_mapping_term()
{
}

template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

Fast_mapping::Fast_mapping(int order)
{
    init(order);
}

template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

std::vector<std::vector<Fast_mapping_terms > > const&
Fast_mapping::get_terms() const
{
    return terms;
}

Fast_mapping::Fast_mapping(Reference_particle const& reference_particle,
        Mapping const& chef_mapping, double mapping_length)
{
    order = chef_mapping.Weight();
    init(order);
    length = mapping_length;
    for (int i = 0; i < 6; ++i) {
        int chef_i = get_chef_index(i);
        Jet__environment_ptr env = chef_mapping.Env();
        for (Jet::const_iterator jet_it = chef_mapping(chef_i).begin(); jet_it
                != chef_mapping(chef_i).end(); ++jet_it) {
            if (jet_it->coefficient() == 0.0) {
                //ignore zero coefficients
            } else {
                int term_order = jet_it->exponents(env).Sum();
                Fast_mapping_term tmp_term(term_order);
                tmp_term.coeff() = jet_it->coefficient();
                int which = 0;
                for (int index = 0; index < 6; ++index) {
                    int chef_i2 = get_chef_index(index);
                    int expt = jet_it->exponents(env)(chef_i2);
                    for (int count = 0; count < expt; ++count) {
                        tmp_term.index(which) = index;
                        ++which;
                    }
                }
                add_term(i, tmp_term);
            }
        }
    }
}

void
Fast_mapping::set_length(double length)
{
    this->length = length;
}

double
Fast_mapping::get_length() const
{
    return length;
}

int
Fast_mapping::get_order() const
{
    return order;
}

void
Fast_mapping::apply(Bunch & bunch)
{

    bunch.get_reference_particle().increment_trajectory(length);
    double temp[6];
    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < local_num; ++part) {
        for (int i = 0; i < 6; ++i) {
            temp[i] = 0.0;
            double term;
            Fast_mapping_terms::const_iterator telem;
            for (telem = terms[i][0].begin(); telem != terms[i][0].end(); ++telem) {
                temp[i] += telem->coeff();
            }
            for (telem = terms[i][1].begin(); telem != terms[i][1].end(); ++telem) {
                term = telem->coeff();
                term *= particles[part][telem->index(0)];
                temp[i] += term;
            }
            if (order > 1) {
                for (telem = terms[i][2].begin(); telem != terms[i][2].end(); ++telem) {
                    term = telem->coeff();
                    term *= particles[part][telem->index(0)];
                    term *= particles[part][telem->index(1)];
                    temp[i] += term;
                }
                if (order > 2) {
                    for (telem = terms[i][3].begin(); telem
                            != terms[i][3].end(); ++telem) {
                        term = telem->coeff();
                        term *= particles[part][telem->index(0)];
                        term *= particles[part][telem->index(1)];
                        term *= particles[part][telem->index(2)];
                        temp[i] += term;
                    }
                    if (order > 3) {
                        for (int suborder = 4; suborder <= order; ++suborder) {
                            for (telem = terms[i][suborder].begin(); telem
                                    != terms[i][suborder].end(); ++telem) {
                                term = telem->coeff();
                                for (int j = 0; j < suborder; ++j) {
                                    term *= particles[part][telem->index(j)];
                                }
                                temp[i] += term;
                            }
                        }
                    }
                }
            }
        }
#ifdef MANUAL_LOOP_UNROLL
        particles[part][0] = temp[0];
        particles[part][1] = temp[1];
        particles[part][2] = temp[2];
        particles[part][3] = temp[3];
        particles[part][4] = temp[4];
        particles[part][5] = temp[5];
#else
        for (int i = 0; i < 6; ++i) {
            particles[ part][i] = temp[i];
        }
#endif
    }

template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
    void
TFast_mapping<double>::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);


template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
