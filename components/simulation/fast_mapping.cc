#include "fast_mapping.h"
#include "components/lattice/chef_utils.h"

#define MANUAL_LOOP_UNROLL shockingly_yes

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

Fast_mapping_term::Fast_mapping_term(int order)
{
    i = new int[order + 1];
    this->the_order = order;
}

Fast_mapping_term::Fast_mapping_term(std::ifstream & stream)
{
    std::stringstream sstream(read_line_ignoring_comments(stream));
    sstream >> the_order;
    i = new int[the_order + 1];
    for (int j = 0; j < the_order + 1; ++j) {
        sstream >> i[j];
    }
    sstream >> the_coeff;
}

Fast_mapping_term::Fast_mapping_term(Fast_mapping_term const& t)
{
    the_coeff = t.the_coeff;
    the_order = t.the_order;
    i = new int[the_order + 1];
    for (int j = 0; j < the_order + 1; ++j) {
        i[j] = t.i[j];
    }
}

void
Fast_mapping_term::write_to_stream(std::ofstream& stream) const
{
    stream << the_order << " ";
    for (int j = 0; j < the_order + 1; ++j) {
        stream << i[j] << " ";
    }
    stream.precision(17);
    stream << the_coeff << std::endl;
}

Fast_mapping_term::~Fast_mapping_term()
{
    delete[] i;
}

void
Fast_mapping::init(int order)
{
    terms.resize(6);
    this->order = order;
    for (int comp_index = 0; comp_index < 6; ++comp_index) {
        terms.at(comp_index).resize(order + 1);
    }
    length = 0.0;
}

Fast_mapping::Fast_mapping(int order)
{
    init(order);
}

Fast_mapping::Fast_mapping(std::string const& filename)
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
                Fast_mapping_term tmp_term(file);
                terms.at(index).at(order).push_back(tmp_term);
            }
        }
    }
    if (read_line_ignoring_comments(file) != "end_fast_mapping") {
        throw std::runtime_error("Fast_mapping::read_from_file(" + filename
                + ") found a truncated file.");
    }
    file.close();
}

inline double
quickpow(double x, int i)
{
    double retval = 1.0;
    while (i--) {
        retval *= x;
    }
    return retval;
}

void
Fast_mapping::add_term(int index, Fast_mapping_term const& term)
{
    terms.at(index).at(term.order()).push_back(term);
}

Fast_mapping::Fast_mapping(Reference_particle const& reference_particle,
        Mapping const& chef_mapping)
{
    std::vector<double > u = chef_unit_conversion(reference_particle);
    order = chef_mapping.Weight();
    init(order);
    for (int i = 0; i < 6; ++i) {
        int chef_i = get_chef_index(i);
        int nterm = 0;
        Jet__environment_ptr env = chef_mapping.Env();
        for (Jet::const_iterator jet_it = chef_mapping(chef_i).begin(); jet_it
                != chef_mapping(chef_i).end(); ++jet_it) {
            if (jet_it->coefficient() == 0.0) {
                //ignore zero coefficients
            } else {
                int term_order = jet_it->exponents(env).Sum();
                Fast_mapping_term tmp_term(term_order);
                tmp_term.coeff() = jet_it->coefficient() * u[i];
                int which = 0;
                for (int index = 0; index < 6; ++index) {
                    int chef_i2 = get_chef_index(index);
                    int expt = jet_it->exponents(env)(chef_i2);
                    for (int count = 0; count < expt; ++count) {
                        tmp_term.index(which) = index;
                        ++which;
                    }
                    tmp_term.coeff() *= 1.0 / quickpow(u[index], expt);
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

void
Fast_mapping::apply(Bunch & bunch)
{
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
}

void
Fast_mapping::write_to_file(std::string const& filename)
{
    std::ofstream file(filename.c_str());
    file << "# Synergia Fast_mapping\n";
    file << "# file format=1.0\n";
    file << "# order:\n";
    file << order << std::endl;
    file << "# length:\n";
    file << length << std::endl;
    for (int index = 0; index < 6; ++index) {
        for (int order = 0; order <= this->order; ++order) {
            file << "# index=" << index << ", order=" << order << std::endl;
            file << "# num:\n";
            file << terms[index][order].size() << std::endl;
            if (terms[index][order].size() > 0) {
                file << "# terms:\n";
                Fast_mapping_terms::const_iterator telem;
                for (telem = terms[index][order].begin(); telem
                        != terms[index][order].end(); ++telem) {
                    telem->write_to_stream(file);
                }
            }
        }
    }
    file << "end_fast_mapping\n";
    file.close();
}
