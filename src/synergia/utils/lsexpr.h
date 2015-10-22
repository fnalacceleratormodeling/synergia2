#ifndef LSEXPR_H
#define LSEXPR_H

#include <string>
#include <vector>
#include <fstream>

struct Lsexpr;

struct Lsexpr
{
    bool has_label;
    std::string label;
    bool is_atom;
    std::string atom;
    std::vector<Lsexpr> list;

    Lsexpr() :
        has_label(false)
      , label("")
      , is_atom(false)
      , atom("")
      , list()
    { }

    Lsexpr(std::string const& atom) :
        has_label(false)
      , label("")
      , is_atom(true)
      , atom(atom)
      , list()
    { }

//    Lsexpr(int atom) :
//        has_label(false)
//      , label("")
//      , is_atom(true)
//      , atom()
//      , list()
//    {
//        std::stringstream ss;
//        ss << atom;
//        this->atom = ss.str();
//    }

//    Lsexpr(double atom) :
//        has_label(false)
//      , label("")
//      , is_atom(true)
//      , atom(atom)
//      , list()
//    {
//        std::stringstream ss;
//        ss << atom;
//        this->atom = ss.str();
//    }

    void set_label(std::string const& label)
    {
        this->label = label;
        has_label = true;
    }
};

Lsexpr read_lsexpr(std::istream & stream);

void write_lsexpr(Lsexpr const& lsexpr, std::ostream & stream);

#endif // LSEXPR_H
