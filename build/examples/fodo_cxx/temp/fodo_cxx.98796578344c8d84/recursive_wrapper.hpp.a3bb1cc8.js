var data = {lines:[
{"lineNum":"    1","line":"//-----------------------------------------------------------------------------"},
{"lineNum":"    2","line":"// boost variant/recursive_wrapper.hpp header file"},
{"lineNum":"    3","line":"// See http://www.boost.org for updates, documentation, and revision history."},
{"lineNum":"    4","line":"//-----------------------------------------------------------------------------"},
{"lineNum":"    5","line":"//"},
{"lineNum":"    6","line":"// Copyright (c) 2002-2003"},
{"lineNum":"    7","line":"// Eric Friedman, Itay Maman"},
{"lineNum":"    8","line":"//"},
{"lineNum":"    9","line":"// Distributed under the Boost Software License, Version 1.0. (See"},
{"lineNum":"   10","line":"// accompanying file LICENSE_1_0.txt or copy at"},
{"lineNum":"   11","line":"// http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"#ifndef BOOST_VARIANT_RECURSIVE_WRAPPER_HPP"},
{"lineNum":"   14","line":"#define BOOST_VARIANT_RECURSIVE_WRAPPER_HPP"},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"#include <boost/variant/recursive_wrapper_fwd.hpp>"},
{"lineNum":"   17","line":"#include <boost/variant/detail/move.hpp>"},
{"lineNum":"   18","line":"#include <boost/checked_delete.hpp>"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"namespace boost {"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"//////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   23","line":"// class template recursive_wrapper"},
{"lineNum":"   24","line":"//"},
{"lineNum":"   25","line":"// See docs and recursive_wrapper_fwd.hpp for more information."},
{"lineNum":"   26","line":"//"},
{"lineNum":"   27","line":""},
{"lineNum":"   28","line":"template <typename T>"},
{"lineNum":"   29","line":"class recursive_wrapper"},
{"lineNum":"   30","line":"{"},
{"lineNum":"   31","line":"public: // typedefs"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"    typedef T type;"},
{"lineNum":"   34","line":""},
{"lineNum":"   35","line":"private: // representation"},
{"lineNum":"   36","line":""},
{"lineNum":"   37","line":"    T* p_;"},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"public: // structors"},
{"lineNum":"   40","line":""},
{"lineNum":"   41","line":"    ~recursive_wrapper();"},
{"lineNum":"   42","line":"    recursive_wrapper();"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"    recursive_wrapper(const recursive_wrapper& operand);"},
{"lineNum":"   45","line":"    recursive_wrapper(const T& operand);"},
{"lineNum":"   46","line":""},
{"lineNum":"   47","line":"#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES"},
{"lineNum":"   48","line":"    recursive_wrapper(recursive_wrapper&& operand);"},
{"lineNum":"   49","line":"    recursive_wrapper(T&& operand);"},
{"lineNum":"   50","line":"#endif"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"private: // helpers, for modifiers (below)"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":"    void assign(const T& rhs);"},
{"lineNum":"   55","line":""},
{"lineNum":"   56","line":"public: // modifiers"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"    recursive_wrapper& operator=(const recursive_wrapper& rhs)"},
{"lineNum":"   59","line":"    {"},
{"lineNum":"   60","line":"        assign( rhs.get() );"},
{"lineNum":"   61","line":"        return *this;"},
{"lineNum":"   62","line":"    }"},
{"lineNum":"   63","line":""},
{"lineNum":"   64","line":"    recursive_wrapper& operator=(const T& rhs)"},
{"lineNum":"   65","line":"    {"},
{"lineNum":"   66","line":"        assign( rhs );"},
{"lineNum":"   67","line":"        return *this;"},
{"lineNum":"   68","line":"    }"},
{"lineNum":"   69","line":""},
{"lineNum":"   70","line":"    void swap(recursive_wrapper& operand) BOOST_NOEXCEPT"},
{"lineNum":"   71","line":"    {"},
{"lineNum":"   72","line":"        T* temp = operand.p_;"},
{"lineNum":"   73","line":"        operand.p_ = p_;"},
{"lineNum":"   74","line":"        p_ = temp;"},
{"lineNum":"   75","line":"    }"},
{"lineNum":"   76","line":""},
{"lineNum":"   77","line":""},
{"lineNum":"   78","line":"#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES"},
{"lineNum":"   79","line":"    recursive_wrapper& operator=(recursive_wrapper&& rhs) BOOST_NOEXCEPT"},
{"lineNum":"   80","line":"    {"},
{"lineNum":"   81","line":"        swap(rhs);"},
{"lineNum":"   82","line":"        return *this;"},
{"lineNum":"   83","line":"    }"},
{"lineNum":"   84","line":""},
{"lineNum":"   85","line":"    recursive_wrapper& operator=(T&& rhs)"},
{"lineNum":"   86","line":"    {"},
{"lineNum":"   87","line":"        get() = detail::variant::move(rhs);"},
{"lineNum":"   88","line":"        return *this;"},
{"lineNum":"   89","line":"    }"},
{"lineNum":"   90","line":"#endif"},
{"lineNum":"   91","line":""},
{"lineNum":"   92","line":"public: // queries"},
{"lineNum":"   93","line":""},
{"lineNum":"   94","line":"    T& get() { return *get_pointer(); }"},
{"lineNum":"   95","line":"    const T& get() const { return *get_pointer(); }"},
{"lineNum":"   96","line":""},
{"lineNum":"   97","line":"    T* get_pointer() { return p_; }","class":"lineNoCov","hits":"0","possible_hits":"15",},
{"lineNum":"   98","line":"    const T* get_pointer() const { return p_; }","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   99","line":""},
{"lineNum":"  100","line":"};"},
{"lineNum":"  101","line":""},
{"lineNum":"  102","line":"template <typename T>"},
{"lineNum":"  103","line":"recursive_wrapper<T>::~recursive_wrapper()"},
{"lineNum":"  104","line":"{"},
{"lineNum":"  105","line":"    boost::checked_delete(p_);","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  106","line":"}"},
{"lineNum":"  107","line":""},
{"lineNum":"  108","line":"template <typename T>"},
{"lineNum":"  109","line":"recursive_wrapper<T>::recursive_wrapper()"},
{"lineNum":"  110","line":"    : p_(new T)"},
{"lineNum":"  111","line":"{"},
{"lineNum":"  112","line":"}"},
{"lineNum":"  113","line":""},
{"lineNum":"  114","line":"template <typename T>"},
{"lineNum":"  115","line":"recursive_wrapper<T>::recursive_wrapper(const recursive_wrapper& operand)"},
{"lineNum":"  116","line":"    : p_(new T( operand.get() ))","class":"lineNoCov","hits":"0","possible_hits":"13",},
{"lineNum":"  117","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  118","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  119","line":""},
{"lineNum":"  120","line":"template <typename T>"},
{"lineNum":"  121","line":"recursive_wrapper<T>::recursive_wrapper(const T& operand)"},
{"lineNum":"  122","line":"    : p_(new T(operand))"},
{"lineNum":"  123","line":"{"},
{"lineNum":"  124","line":"}"},
{"lineNum":"  125","line":""},
{"lineNum":"  126","line":"#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES"},
{"lineNum":"  127","line":"template <typename T>"},
{"lineNum":"  128","line":"recursive_wrapper<T>::recursive_wrapper(recursive_wrapper&& operand)"},
{"lineNum":"  129","line":"    : p_(new T( detail::variant::move(operand.get()) ))","class":"lineNoCov","hits":"0","possible_hits":"11",},
{"lineNum":"  130","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  131","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  132","line":""},
{"lineNum":"  133","line":"template <typename T>"},
{"lineNum":"  134","line":"recursive_wrapper<T>::recursive_wrapper(T&& operand)"},
{"lineNum":"  135","line":"    : p_(new T( detail::variant::move(operand) ))","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  136","line":"{"},
{"lineNum":"  137","line":"}"},
{"lineNum":"  138","line":"#endif"},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":"template <typename T>"},
{"lineNum":"  141","line":"void recursive_wrapper<T>::assign(const T& rhs)"},
{"lineNum":"  142","line":"{"},
{"lineNum":"  143","line":"    this->get() = rhs;"},
{"lineNum":"  144","line":"}"},
{"lineNum":"  145","line":""},
{"lineNum":"  146","line":"// function template swap"},
{"lineNum":"  147","line":"//"},
{"lineNum":"  148","line":"// Swaps two recursive_wrapper<T> objects of the same type T."},
{"lineNum":"  149","line":"//"},
{"lineNum":"  150","line":"template <typename T>"},
{"lineNum":"  151","line":"inline void swap(recursive_wrapper<T>& lhs, recursive_wrapper<T>& rhs) BOOST_NOEXCEPT"},
{"lineNum":"  152","line":"{"},
{"lineNum":"  153","line":"    lhs.swap(rhs);"},
{"lineNum":"  154","line":"}"},
{"lineNum":"  155","line":""},
{"lineNum":"  156","line":"} // namespace boost"},
{"lineNum":"  157","line":""},
{"lineNum":"  158","line":"#endif // BOOST_VARIANT_RECURSIVE_WRAPPER_HPP"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 10, "covered" : 0,};
var merged_data = [];
