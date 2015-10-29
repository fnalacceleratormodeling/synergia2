#ifndef LSEXPR_H
#define LSEXPR_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

struct Lsexpr;

struct Lsexpr
{
    typedef std::vector<Lsexpr> sequence_t;
    typedef sequence_t::const_iterator const_iterator_t;
    typedef sequence_t::iterator iterator_t;

    bool has_label;
    std::string label;
    bool is_atom;
    std::string atom;
    sequence_t sequence;

    Lsexpr() :
        has_label(false)
      , label("")
      , is_atom(false)
      , atom("")
      , sequence()
    { }

    Lsexpr(std::string const& atom) :
        has_label(false)
      , label("")
      , is_atom(true)
      , atom(atom)
      , sequence()
    { }

    std::string to_string(int i) const
    {
        std::stringstream ss;
        ss << i;
        return ss.str();
    }

    int stoi(std::string const& s) const
    {
        int retval;
        std::istringstream ss(s);
        ss >> retval;
        return retval;
    }

    std::string to_string(double d) const
    {
        std::stringstream ss;
        const int precision = 15;
        ss.precision(precision);
        ss << d;
        return(ss.str());
    }

    double stod(std::string const& s) const
    {
        double retval;
        std::istringstream ss(s);
        ss >> retval;
        return retval;
    }

    Lsexpr(int atom) :
        has_label(false)
      , label("")
      , is_atom(true)
      , atom(to_string(atom))
      , sequence()
    {
    }

    Lsexpr(double atom) :
        has_label(false)
      , label("")
      , is_atom(true)
      , atom(to_string(atom))
      , sequence()
    {
    }

    Lsexpr(std::vector<std::string> const& sequence) :
        has_label(false)
      , label("")
      , is_atom(false)
      , atom()
      , sequence()
    {
        for(std::vector<std::string>::const_iterator it = sequence.begin();
            it != sequence.end(); ++it) {
            this->sequence.push_back(Lsexpr(*it));
        }
    }

    Lsexpr(std::vector<int> const& sequence) :
        has_label(false)
      , label("")
      , is_atom(false)
      , atom()
      , sequence()
    {
        for(std::vector<int>::const_iterator it = sequence.begin();
            it != sequence.end(); ++it) {
            this->sequence.push_back(Lsexpr(*it));
        }
    }

    Lsexpr(std::vector<double> const& sequence) :
        has_label(false)
      , label("")
      , is_atom(false)
      , atom()
      , sequence()
    {
        for(std::vector<double>::const_iterator it = sequence.begin();
            it != sequence.end(); ++it) {
            this->sequence.push_back(Lsexpr(*it));
        }
    }

    void set_label(std::string const& label)
    {
        this->label = label;
        has_label = true;
    }

    bool is_labeled() const
    {
        return has_label;
    }

    std::string const& get_label() const
    {
        return label;
    }

    bool is_atomic() const
    {
        return is_atom;
    }

    std::string get_string() const
    {
        if(!is_atom) {
            throw std::runtime_error("Lsexpr::get_string: non-atomic Lsexpr");
        }
        return atom;
    }

    int get_int() const
    {
        if(!is_atom) {
            throw std::runtime_error("Lsexpr::get_int: non-atomic Lsexpr");
        }
        return stoi(atom);
    }

    double get_double() const
    {
        if(!is_atom) {
            throw std::runtime_error("Lsexpr::get_double: non-atomic Lsexpr");
        }
        return stod(atom);
    }

    void push_back(Lsexpr const& lsexpr)
    {
        sequence.push_back(lsexpr);
    }

    const_iterator_t begin() const
    {
        return sequence.begin();
    }

    const_iterator_t end() const
    {
        return sequence.end();
    }

    iterator_t begin()
    {
        return sequence.begin();
    }

    iterator_t end()
    {
        return sequence.end();
    }

    std::vector<std::string > get_string_vector() const
    {
        if(is_atom) {
            throw std::runtime_error("Lsexpr::get_string_vector: atomic Lsexpr");
        }
        std::vector<std::string> retval;
        for(const_iterator_t it = begin(); it!= end();
            ++it) {
            retval.push_back(it->get_string());
        }
        return retval;
    }

    std::vector<int> get_int_vector() const
    {
        if(is_atom) {
            throw std::runtime_error("Lsexpr::get_int_vector: atomic Lsexpr");
        }
        std::vector<int> retval;
        for(const_iterator_t it = begin(); it!= end();
            ++it) {
            retval.push_back(stoi(it->get_string()));
        }
        return retval;
    }

    std::vector<double> get_double_vector() const
    {
        if(is_atom) {
            throw std::runtime_error("Lsexpr::get_double_vector: atomic Lsexpr");
        }
        std::vector<double> retval;
        for(const_iterator_t it = begin(); it!= end();
            ++it) {
            retval.push_back(stod(it->get_string()));
        }
        return retval;
    }

};

Lsexpr read_lsexpr(std::istream & stream);

void write_lsexpr(Lsexpr const& lsexpr, std::ostream & stream);

#endif // LSEXPR_H
