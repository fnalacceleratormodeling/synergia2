
#include "synergia/lattice/chef_utils.h"
#define MANUAL_LOOP_UNROLL shockingly_yes

std::string
read_line_ignoring_comments(std::ifstream &file);

template <typename T>
Fast_mapping_term<T>::Fast_mapping_term(int order) : i(order + 1)
{
    this->the_order = order;
}

template <typename T>
Fast_mapping_term<T>::Fast_mapping_term(std::ifstream & stream)
{
    std::stringstream sstream(read_line_ignoring_comments(stream));
    sstream >> the_order;
    i.resize(the_order + 1);
    for (int j = 0; j < the_order + 1; ++j) {
        sstream >> i[j];
    }
    sstream >> the_coeff;
}

template <typename T>
Fast_mapping_term<T>::Fast_mapping_term()
{
}

template <typename T>
Fast_mapping_term<T>::Fast_mapping_term(Fast_mapping_term const& t) :
        i(t.the_order + 1)
{
    the_coeff = t.the_coeff;
    the_order = t.the_order;
    for (int j = 0; j < the_order + 1; ++j) {
        i[j] = t.i[j];
    }
}

template <typename T>
void
Fast_mapping_term<T>::write_to_stream(std::ostream& stream) const
{
    stream << the_order << " ";
    for (int j = 0; j < the_order + 1; ++j) {
        stream << i[j] << " ";
    }
    stream.precision(17);
    stream << the_coeff << std::endl;
}

template<typename T>
template<typename Archive>
    void
    Fast_mapping_term<T>::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(the_coeff);
        ar & BOOST_SERIALIZATION_NVP(the_order);
        ar & BOOST_SERIALIZATION_NVP(i);
    }



template<typename T>
Fast_mapping_term<T>::~Fast_mapping_term()
{
}


template<typename T>
void
TFast_mapping<T>::init(int order)
{
    terms.resize(6);
    this->order = order;
    for (int comp_index = 0; comp_index < 6; ++comp_index) {
        terms.at(comp_index).resize(order + 1);
    }
    length = 0.0;
}

template<typename T>
TFast_mapping<T>::TFast_mapping(int order)
{
    init(order);
}

template<typename T>
TFast_mapping<T>::TFast_mapping(std::string const& filename)
{
    std::ifstream file(filename.c_str());
    std::stringstream stream(read_line_ignoring_comments(file));
    stream >> order;
    init(order);
    stream.str(read_line_ignoring_comments(file));
    stream >> length;
    int num_terms_read;
    for (int index = 0; index < 6; ++index) {
        for (int order = 0; order <= this->order; ++order) {
            stream.str(read_line_ignoring_comments(file));
            stream.clear();
            stream >> num_terms_read;
            //            std::cout << "jfa: num_terms_read = " << num_terms_read << std::endl;
            for (int term_index = 0; term_index < num_terms_read; ++term_index) {
                Fast_mapping_term<T> tmp_term(file);
                terms.at(index).at(order).push_back(tmp_term);
            }
        }
    }
    if (read_line_ignoring_comments(file) != "end_fast_mapping") {
        throw std::runtime_error(
                "Fast_mapping::read_from_file(" + filename
                        + ") found a truncated file.");
    }
    file.close();
}

template<typename T>
TFast_mapping<T>::TFast_mapping()
{
}

template<typename T>
inline T
quickpow(T x, int i)
{
    T retval = 1.0;
    while (i--) {
        retval *= x;
    }
    return retval;
}

template<typename T>
void
TFast_mapping<T>::add_term(int index, Fast_mapping_term<T> const& term)
{
    terms.at(index).at(term.order()).push_back(term);
}

template<typename T>
std::vector<std::vector<std::list<Fast_mapping_term<T> > > >const&
TFast_mapping<T>::get_terms() const
{
    return terms;
}

template<typename T>
TFast_mapping<T>::TFast_mapping(TMapping<T> const& chef_mapping, double mapping_length)

{        
    order = chef_mapping.Weight();  
    init(order);   
    length = mapping_length;
    for (int i = 0; i < 6; ++i) {
        int chef_i = get_chef_index(i);
        EnvPtr<T> env = chef_mapping.Env();
        for (typename TJet<T>::const_iterator jet_it = chef_mapping(chef_i).begin(); jet_it
                != chef_mapping(chef_i).end(); ++jet_it) {
             if (jet_it->coefficient() == T()) {
                 //ignore zero coefficients             
             } else {
                  
                  int term_order = jet_it->exponents(env).Sum();    
                  Fast_mapping_term<T> tmp_term(term_order);
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
                  if (i==0){
              
                  }
                  add_term(i, tmp_term);
              }                
         }
    } 
}

template<typename T>
void
TFast_mapping<T>::set_length(double length)
{
    this->length = length;
}

template<typename T>
double
TFast_mapping<T>::get_length() const
{
    return length;
}

template<typename T>
int
TFast_mapping<T>::get_order() const
{
    return order;
}

template<typename T>
void
TFast_mapping<T>::apply(Bunch & bunch)
{

    bunch.get_reference_particle().increment_trajectory(length);
    double temp[6];
    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < local_num; ++part) {
        for (int i = 0; i < 6; ++i) {
            temp[i] = 0.0;
            double term;
            typename std::list<Fast_mapping_term<T> >::const_iterator telem;
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


}

template<typename T>
void
TFast_mapping<T>::apply(boost::multi_array_ref<T, 1 > coords){
  
   T temp[6];
   for (int i = 0; i < 6; ++i) {
            temp[i] = T();
            T term;
            typename std::list<Fast_mapping_term<T> >::const_iterator telem;
            for (telem = terms[i][0].begin(); telem != terms[i][0].end(); ++telem) {
                temp[i] += telem->coeff();
            }
            for (telem = terms[i][1].begin(); telem != terms[i][1].end(); ++telem) {
                term = telem->coeff();
                term *= coords[telem->index(0)];
                temp[i] += term;
            }
             if (order > 1) {
                for (telem = terms[i][2].begin(); telem != terms[i][2].end(); ++telem) {
                    term = telem->coeff();
                    term *= coords[telem->index(0)];
                    term *= coords[telem->index(1)];
                    temp[i] += term;
                }
                 if (order > 2) {
                    for (telem = terms[i][3].begin(); telem
                            != terms[i][3].end(); ++telem) {
                        term = telem->coeff();
                        term *= coords[telem->index(0)];
                        term *= coords[telem->index(1)];
                        term *= coords[telem->index(2)];
                        temp[i] += term;
                    }
                    if (order > 3) {
                        for (int suborder = 4; suborder <= order; ++suborder) {
                            for (telem = terms[i][suborder].begin(); telem
                                    != terms[i][suborder].end(); ++telem) {
                                term = telem->coeff();
                                for (int j = 0; j < suborder; ++j) {
                                    term *= coords[telem->index(j)];
                                }
                                temp[i] += term;
                            }
                        }
                    }
                 }
             }
    }
    for (int i = 0; i < 6; ++i) {
           coords[i] = temp[i];
    }
  
}


template<typename T>
std::string
TFast_mapping<T>::as_string() const
{
    std::stringstream sstream;
    sstream << "# order:\n";
    sstream << order << std::endl;
    sstream << "# length:\n";
    sstream << length << std::endl;
    for (int index = 0; index < 6; ++index) {
        for (int order = 0; order <= this->order; ++order) {
            sstream << "# index=" << index << ", order=" << order << std::endl;
            sstream << "# num:\n";
            sstream << terms[index][order].size() << std::endl;
            if (terms[index][order].size() > 0) {
                sstream << "# terms:\n";
                typename std::list<Fast_mapping_term<T> > ::const_iterator telem;
                for (telem = terms[index][order].begin(); telem
                        != terms[index][order].end(); ++telem) {
                    telem->write_to_stream(sstream);
                }
            }
        }
    }
    return sstream.str();
}

template<typename T>
void
TFast_mapping<T>::write_to_file(std::string const& filename)
{
    std::ofstream file(filename.c_str());
    file << "# Synergia Fast_mapping\n";
    file << "# file format=1.0\n";
    file << as_string();
    file << "end_fast_mapping\n";
    file.close();
}



template<typename T>
template<typename Archive>
    void
    TFast_mapping<T>::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(terms);
        ar & BOOST_SERIALIZATION_NVP(order);
        ar & BOOST_SERIALIZATION_NVP(length);
    }


