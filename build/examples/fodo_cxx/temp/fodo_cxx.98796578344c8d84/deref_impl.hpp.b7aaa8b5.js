var data = {lines:[
{"lineNum":"    1","line":"/*============================================================================="},
{"lineNum":"    2","line":"    Copyright (c) 2001-2011 Joel de Guzman"},
{"lineNum":"    3","line":"    Copyright (c) 2005 Eric Niebler"},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":"    Distributed under the Boost Software License, Version 1.0. (See accompanying"},
{"lineNum":"    6","line":"    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    7","line":"==============================================================================*/"},
{"lineNum":"    8","line":"#if !defined(FUSION_DEREF_IMPL_07172005_0831)"},
{"lineNum":"    9","line":"#define FUSION_DEREF_IMPL_07172005_0831"},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"#include <boost/fusion/support/config.hpp>"},
{"lineNum":"   12","line":"#include <boost/mpl/eval_if.hpp>"},
{"lineNum":"   13","line":"#include <boost/type_traits/is_const.hpp>"},
{"lineNum":"   14","line":"#include <boost/type_traits/add_const.hpp>"},
{"lineNum":"   15","line":"#include <boost/type_traits/add_reference.hpp>"},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"namespace boost { namespace fusion"},
{"lineNum":"   18","line":"{"},
{"lineNum":"   19","line":"    struct cons_iterator_tag;"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"    namespace extension"},
{"lineNum":"   22","line":"    {"},
{"lineNum":"   23","line":"        template <typename Tag>"},
{"lineNum":"   24","line":"        struct deref_impl;"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"        template <>"},
{"lineNum":"   27","line":"        struct deref_impl<cons_iterator_tag>"},
{"lineNum":"   28","line":"        {"},
{"lineNum":"   29","line":"            template <typename Iterator>"},
{"lineNum":"   30","line":"            struct apply"},
{"lineNum":"   31","line":"            {"},
{"lineNum":"   32","line":"                typedef typename Iterator::cons_type cons_type;"},
{"lineNum":"   33","line":"                typedef typename cons_type::car_type value_type;"},
{"lineNum":"   34","line":""},
{"lineNum":"   35","line":"                typedef typename mpl::eval_if<"},
{"lineNum":"   36","line":"                    is_const<cons_type>"},
{"lineNum":"   37","line":"                  , add_reference<typename add_const<value_type>::type>"},
{"lineNum":"   38","line":"                  , add_reference<value_type> >::type"},
{"lineNum":"   39","line":"                type;"},
{"lineNum":"   40","line":""},
{"lineNum":"   41","line":"                BOOST_CONSTEXPR BOOST_FUSION_GPU_ENABLED"},
{"lineNum":"   42","line":"                static type"},
{"lineNum":"   43","line":"                call(Iterator const& i)"},
{"lineNum":"   44","line":"                {"},
{"lineNum":"   45","line":"                    return i.cons.car;","class":"lineNoCov","hits":"0","possible_hits":"66",},
{"lineNum":"   46","line":"                }"},
{"lineNum":"   47","line":"            };"},
{"lineNum":"   48","line":"        };"},
{"lineNum":"   49","line":"    }"},
{"lineNum":"   50","line":"}}"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"#endif"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":""},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 1, "covered" : 0,};
var merged_data = [];
