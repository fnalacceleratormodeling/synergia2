var data = {lines:[
{"lineNum":"    1","line":"/*============================================================================="},
{"lineNum":"    2","line":"    Copyright (c) 2001-2011 Joel de Guzman"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":"    Distributed under the Boost Software License, Version 1.0. (See accompanying"},
{"lineNum":"    5","line":"    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    6","line":"=============================================================================*/"},
{"lineNum":"    7","line":"#if !defined(SPIRIT_NOT_PREDICATE_MARCH_23_2007_0618PM)"},
{"lineNum":"    8","line":"#define SPIRIT_NOT_PREDICATE_MARCH_23_2007_0618PM"},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#if defined(_MSC_VER)"},
{"lineNum":"   11","line":"#pragma once"},
{"lineNum":"   12","line":"#endif"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"#include <boost/spirit/home/qi/domain.hpp>"},
{"lineNum":"   15","line":"#include <boost/spirit/home/qi/meta_compiler.hpp>"},
{"lineNum":"   16","line":"#include <boost/spirit/home/qi/parser.hpp>"},
{"lineNum":"   17","line":"#include <boost/spirit/home/qi/detail/attributes.hpp>"},
{"lineNum":"   18","line":"#include <boost/spirit/home/support/has_semantic_action.hpp>"},
{"lineNum":"   19","line":"#include <boost/spirit/home/support/handles_container.hpp>"},
{"lineNum":"   20","line":"#include <boost/spirit/home/support/info.hpp>"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"namespace boost { namespace spirit"},
{"lineNum":"   23","line":"{"},
{"lineNum":"   24","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   25","line":"    // Enablers"},
{"lineNum":"   26","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   27","line":"    template <>"},
{"lineNum":"   28","line":"    struct use_operator<qi::domain, proto::tag::logical_not> // enables !p"},
{"lineNum":"   29","line":"      : mpl::true_ {};"},
{"lineNum":"   30","line":"}}"},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":"namespace boost { namespace spirit { namespace qi"},
{"lineNum":"   33","line":"{"},
{"lineNum":"   34","line":"    template <typename Subject>"},
{"lineNum":"   35","line":"    struct not_predicate : unary_parser<not_predicate<Subject> >"},
{"lineNum":"   36","line":"    {"},
{"lineNum":"   37","line":"        typedef Subject subject_type;"},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"        template <typename Context, typename Iterator>"},
{"lineNum":"   40","line":"        struct attribute"},
{"lineNum":"   41","line":"        {"},
{"lineNum":"   42","line":"            typedef unused_type type;"},
{"lineNum":"   43","line":"        };"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"        not_predicate(Subject const& subject_)"},
{"lineNum":"   46","line":"          : subject(subject_) {}"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"        template <typename Iterator, typename Context"},
{"lineNum":"   49","line":"          , typename Skipper, typename Attribute>"},
{"lineNum":"   50","line":"        bool parse(Iterator& first, Iterator const& last"},
{"lineNum":"   51","line":"          , Context& context, Skipper const& skipper"},
{"lineNum":"   52","line":"          , Attribute& /*attr*/) const"},
{"lineNum":"   53","line":"        {"},
{"lineNum":"   54","line":"            Iterator i = first;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   55","line":"            return !subject.parse(i, last, context, skipper, unused);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   56","line":"        }"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"        template <typename Context>"},
{"lineNum":"   59","line":"        info what(Context& context) const"},
{"lineNum":"   60","line":"        {"},
{"lineNum":"   61","line":"            return info(\"not-predicate\", subject.what(context));"},
{"lineNum":"   62","line":"        }"},
{"lineNum":"   63","line":""},
{"lineNum":"   64","line":"        Subject subject;"},
{"lineNum":"   65","line":"    };"},
{"lineNum":"   66","line":""},
{"lineNum":"   67","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   68","line":"    // Parser generators: make_xxx function (objects)"},
{"lineNum":"   69","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   70","line":"    template <typename Elements, typename Modifiers>"},
{"lineNum":"   71","line":"    struct make_composite<proto::tag::logical_not, Elements, Modifiers>"},
{"lineNum":"   72","line":"      : make_unary_composite<Elements, not_predicate>"},
{"lineNum":"   73","line":"    {};"},
{"lineNum":"   74","line":"}}}"},
{"lineNum":"   75","line":""},
{"lineNum":"   76","line":"namespace boost { namespace spirit { namespace traits"},
{"lineNum":"   77","line":"{"},
{"lineNum":"   78","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   79","line":"    template <typename Subject>"},
{"lineNum":"   80","line":"    struct has_semantic_action<qi::not_predicate<Subject> >"},
{"lineNum":"   81","line":"      : unary_has_semantic_action<Subject> {};"},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   84","line":"    template <typename Subject, typename Attribute, typename Context"},
{"lineNum":"   85","line":"        , typename Iterator>"},
{"lineNum":"   86","line":"    struct handles_container<qi::not_predicate<Subject>, Attribute"},
{"lineNum":"   87","line":"        , Context, Iterator>"},
{"lineNum":"   88","line":"      : unary_handles_container<Subject, Attribute, Context, Iterator> {};"},
{"lineNum":"   89","line":"}}}"},
{"lineNum":"   90","line":""},
{"lineNum":"   91","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 2, "covered" : 0,};
var merged_data = [];
