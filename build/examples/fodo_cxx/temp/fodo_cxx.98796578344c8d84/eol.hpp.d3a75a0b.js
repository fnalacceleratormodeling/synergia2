var data = {lines:[
{"lineNum":"    1","line":"/*============================================================================="},
{"lineNum":"    2","line":"    Copyright (c) 2001-2011 Hartmut Kaiser"},
{"lineNum":"    3","line":"    Copyright (c) 2001-2011 Joel de Guzman"},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":"    Distributed under the Boost Software License, Version 1.0. (See accompanying"},
{"lineNum":"    6","line":"    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    7","line":"==============================================================================*/"},
{"lineNum":"    8","line":"#if !defined(BOOST_SPIRIT_EOL_APRIL_18_2008_0751PM)"},
{"lineNum":"    9","line":"#define BOOST_SPIRIT_EOL_APRIL_18_2008_0751PM"},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"#if defined(_MSC_VER)"},
{"lineNum":"   12","line":"#pragma once"},
{"lineNum":"   13","line":"#endif"},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"#include <boost/mpl/bool.hpp>"},
{"lineNum":"   16","line":"#include <boost/spirit/home/qi/domain.hpp>"},
{"lineNum":"   17","line":"#include <boost/spirit/home/qi/parser.hpp>"},
{"lineNum":"   18","line":"#include <boost/spirit/home/qi/meta_compiler.hpp>"},
{"lineNum":"   19","line":"#include <boost/spirit/home/qi/skip_over.hpp>"},
{"lineNum":"   20","line":"#include <boost/spirit/home/support/common_terminals.hpp>"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"namespace boost { namespace spirit"},
{"lineNum":"   23","line":"{"},
{"lineNum":"   24","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   25","line":"    // Enablers"},
{"lineNum":"   26","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   27","line":"    template <>"},
{"lineNum":"   28","line":"    struct use_terminal<qi::domain, tag::eol>       // enables eol"},
{"lineNum":"   29","line":"      : mpl::true_ {};"},
{"lineNum":"   30","line":"}}"},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":"namespace boost { namespace spirit { namespace qi"},
{"lineNum":"   33","line":"{"},
{"lineNum":"   34","line":"#ifndef BOOST_SPIRIT_NO_PREDEFINED_TERMINALS"},
{"lineNum":"   35","line":"    using spirit::eol;"},
{"lineNum":"   36","line":"#endif"},
{"lineNum":"   37","line":"    using spirit::eol_type;"},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"    struct eol_parser : primitive_parser<eol_parser>"},
{"lineNum":"   40","line":"    {"},
{"lineNum":"   41","line":"        template <typename Context, typename Iterator>"},
{"lineNum":"   42","line":"        struct attribute"},
{"lineNum":"   43","line":"        {"},
{"lineNum":"   44","line":"            typedef unused_type type;"},
{"lineNum":"   45","line":"        };"},
{"lineNum":"   46","line":""},
{"lineNum":"   47","line":"        template <typename Iterator, typename Context"},
{"lineNum":"   48","line":"          , typename Skipper, typename Attribute>"},
{"lineNum":"   49","line":"        bool parse(Iterator& first, Iterator const& last"},
{"lineNum":"   50","line":"          , Context& /*context*/, Skipper const& skipper"},
{"lineNum":"   51","line":"          , Attribute& /*attr*/) const"},
{"lineNum":"   52","line":"        {"},
{"lineNum":"   53","line":"            qi::skip_over(first, last, skipper);"},
{"lineNum":"   54","line":""},
{"lineNum":"   55","line":"            Iterator it = first;"},
{"lineNum":"   56","line":"            bool matched = false;"},
{"lineNum":"   57","line":"            if (it != last && *it == \'\\r\')  // CR","class":"lineNoCov","hits":"0","possible_hits":"21",},
{"lineNum":"   58","line":"            {"},
{"lineNum":"   59","line":"                matched = true;"},
{"lineNum":"   60","line":"                ++it;"},
{"lineNum":"   61","line":"            }"},
{"lineNum":"   62","line":"            if (it != last && *it == \'\\n\')  // LF","class":"lineNoCov","hits":"0","possible_hits":"10",},
{"lineNum":"   63","line":"            {"},
{"lineNum":"   64","line":"                matched = true;"},
{"lineNum":"   65","line":"                ++it;"},
{"lineNum":"   66","line":"            }"},
{"lineNum":"   67","line":""},
{"lineNum":"   68","line":"            if (!matched)","class":"lineNoCov","hits":"0","possible_hits":"10",},
{"lineNum":"   69","line":"                return false;"},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"            first = it;"},
{"lineNum":"   72","line":"            return true;"},
{"lineNum":"   73","line":"        }"},
{"lineNum":"   74","line":""},
{"lineNum":"   75","line":"        template <typename Context>"},
{"lineNum":"   76","line":"        info what(Context& /*context*/) const"},
{"lineNum":"   77","line":"        {"},
{"lineNum":"   78","line":"            return info(\"eol\");"},
{"lineNum":"   79","line":"        }"},
{"lineNum":"   80","line":"    };"},
{"lineNum":"   81","line":""},
{"lineNum":"   82","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   83","line":"    // Parser generators: make_xxx function (objects)"},
{"lineNum":"   84","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   85","line":"    template <typename Modifiers>"},
{"lineNum":"   86","line":"    struct make_primitive<tag::eol, Modifiers>"},
{"lineNum":"   87","line":"    {"},
{"lineNum":"   88","line":"        typedef eol_parser result_type;"},
{"lineNum":"   89","line":"        result_type operator()(unused_type, unused_type) const"},
{"lineNum":"   90","line":"        {"},
{"lineNum":"   91","line":"            return result_type();"},
{"lineNum":"   92","line":"        }"},
{"lineNum":"   93","line":"    };"},
{"lineNum":"   94","line":"}}}"},
{"lineNum":"   95","line":""},
{"lineNum":"   96","line":"#endif"},
{"lineNum":"   97","line":""},
{"lineNum":"   98","line":""},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 3, "covered" : 0,};
var merged_data = [];
