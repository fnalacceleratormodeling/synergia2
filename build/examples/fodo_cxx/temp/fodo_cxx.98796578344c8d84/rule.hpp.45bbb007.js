var data = {lines:[
{"lineNum":"    1","line":"/*============================================================================="},
{"lineNum":"    2","line":"    Copyright (c) 2001-2011 Joel de Guzman"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":"    Distributed under the Boost Software License, Version 1.0. (See accompanying"},
{"lineNum":"    5","line":"    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    6","line":"==============================================================================*/"},
{"lineNum":"    7","line":"#if !defined(BOOST_SPIRIT_RULE_FEBRUARY_12_2007_1020AM)"},
{"lineNum":"    8","line":"#define BOOST_SPIRIT_RULE_FEBRUARY_12_2007_1020AM"},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#if defined(_MSC_VER)"},
{"lineNum":"   11","line":"#pragma once"},
{"lineNum":"   12","line":"#endif"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"#include <boost/assert.hpp>"},
{"lineNum":"   15","line":"#include <boost/static_assert.hpp>"},
{"lineNum":"   16","line":"#include <boost/config.hpp>"},
{"lineNum":"   17","line":"#include <boost/function.hpp>"},
{"lineNum":"   18","line":"#include <boost/mpl/vector.hpp>"},
{"lineNum":"   19","line":"#include <boost/type_traits/add_reference.hpp>"},
{"lineNum":"   20","line":"#include <boost/type_traits/is_convertible.hpp>"},
{"lineNum":"   21","line":"#include <boost/type_traits/is_same.hpp>"},
{"lineNum":"   22","line":""},
{"lineNum":"   23","line":"#include <boost/fusion/include/vector.hpp>"},
{"lineNum":"   24","line":"#include <boost/fusion/include/size.hpp>"},
{"lineNum":"   25","line":"#include <boost/fusion/include/make_vector.hpp>"},
{"lineNum":"   26","line":"#include <boost/fusion/include/cons.hpp>"},
{"lineNum":"   27","line":"#include <boost/fusion/include/as_list.hpp>"},
{"lineNum":"   28","line":"#include <boost/fusion/include/as_vector.hpp>"},
{"lineNum":"   29","line":""},
{"lineNum":"   30","line":"#include <boost/spirit/home/support/unused.hpp>"},
{"lineNum":"   31","line":"#include <boost/spirit/home/support/argument.hpp>"},
{"lineNum":"   32","line":"#include <boost/spirit/home/support/context.hpp>"},
{"lineNum":"   33","line":"#include <boost/spirit/home/support/info.hpp>"},
{"lineNum":"   34","line":"#include <boost/spirit/home/qi/detail/attributes.hpp>"},
{"lineNum":"   35","line":"#include <boost/spirit/home/support/nonterminal/extract_param.hpp>"},
{"lineNum":"   36","line":"#include <boost/spirit/home/support/nonterminal/locals.hpp>"},
{"lineNum":"   37","line":"#include <boost/spirit/home/qi/reference.hpp>"},
{"lineNum":"   38","line":"#include <boost/spirit/home/qi/nonterminal/detail/parameterized.hpp>"},
{"lineNum":"   39","line":"#include <boost/spirit/home/qi/nonterminal/detail/parser_binder.hpp>"},
{"lineNum":"   40","line":"#include <boost/spirit/home/qi/nonterminal/nonterminal_fwd.hpp>"},
{"lineNum":"   41","line":"#include <boost/spirit/home/qi/skip_over.hpp>"},
{"lineNum":"   42","line":""},
{"lineNum":"   43","line":"#if defined(BOOST_MSVC)"},
{"lineNum":"   44","line":"# pragma warning(push)"},
{"lineNum":"   45","line":"# pragma warning(disable: 4355) // \'this\' : used in base member initializer list warning"},
{"lineNum":"   46","line":"# pragma warning(disable: 4127) // conditional expression is constant"},
{"lineNum":"   47","line":"#endif"},
{"lineNum":"   48","line":""},
{"lineNum":"   49","line":"namespace boost { namespace spirit { namespace qi"},
{"lineNum":"   50","line":"{"},
{"lineNum":"   51","line":"    BOOST_PP_REPEAT(SPIRIT_ATTRIBUTES_LIMIT, SPIRIT_USING_ATTRIBUTE, _)"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"    using spirit::_pass_type;"},
{"lineNum":"   54","line":"    using spirit::_val_type;"},
{"lineNum":"   55","line":"    using spirit::_a_type;"},
{"lineNum":"   56","line":"    using spirit::_b_type;"},
{"lineNum":"   57","line":"    using spirit::_c_type;"},
{"lineNum":"   58","line":"    using spirit::_d_type;"},
{"lineNum":"   59","line":"    using spirit::_e_type;"},
{"lineNum":"   60","line":"    using spirit::_f_type;"},
{"lineNum":"   61","line":"    using spirit::_g_type;"},
{"lineNum":"   62","line":"    using spirit::_h_type;"},
{"lineNum":"   63","line":"    using spirit::_i_type;"},
{"lineNum":"   64","line":"    using spirit::_j_type;"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"#ifndef BOOST_SPIRIT_NO_PREDEFINED_TERMINALS"},
{"lineNum":"   67","line":""},
{"lineNum":"   68","line":"    using spirit::_pass;"},
{"lineNum":"   69","line":"    using spirit::_val;"},
{"lineNum":"   70","line":"    using spirit::_a;"},
{"lineNum":"   71","line":"    using spirit::_b;"},
{"lineNum":"   72","line":"    using spirit::_c;"},
{"lineNum":"   73","line":"    using spirit::_d;"},
{"lineNum":"   74","line":"    using spirit::_e;"},
{"lineNum":"   75","line":"    using spirit::_f;"},
{"lineNum":"   76","line":"    using spirit::_g;"},
{"lineNum":"   77","line":"    using spirit::_h;"},
{"lineNum":"   78","line":"    using spirit::_i;"},
{"lineNum":"   79","line":"    using spirit::_j;"},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"#endif"},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"    using spirit::info;"},
{"lineNum":"   84","line":"    using spirit::locals;"},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"    template <"},
{"lineNum":"   87","line":"        typename Iterator, typename T1, typename T2, typename T3"},
{"lineNum":"   88","line":"      , typename T4>"},
{"lineNum":"   89","line":"    struct rule","class":"lineNoCov","hits":"0","possible_hits":"51",},
{"lineNum":"   90","line":"      : proto::extends<"},
{"lineNum":"   91","line":"            typename proto::terminal<"},
{"lineNum":"   92","line":"                reference<rule<Iterator, T1, T2, T3, T4> const>"},
{"lineNum":"   93","line":"            >::type"},
{"lineNum":"   94","line":"          , rule<Iterator, T1, T2, T3, T4>"},
{"lineNum":"   95","line":"        >"},
{"lineNum":"   96","line":"      , parser<rule<Iterator, T1, T2, T3, T4> >"},
{"lineNum":"   97","line":"    {"},
{"lineNum":"   98","line":"        typedef Iterator iterator_type;"},
{"lineNum":"   99","line":"        typedef rule<Iterator, T1, T2, T3, T4> this_type;"},
{"lineNum":"  100","line":"        typedef reference<this_type const> reference_;"},
{"lineNum":"  101","line":"        typedef typename proto::terminal<reference_>::type terminal;"},
{"lineNum":"  102","line":"        typedef proto::extends<terminal, this_type> base_type;"},
{"lineNum":"  103","line":"        typedef mpl::vector<T1, T2, T3, T4> template_params;"},
{"lineNum":"  104","line":""},
{"lineNum":"  105","line":"        // The rule\'s locals_type: a sequence of types to be used as local variables"},
{"lineNum":"  106","line":"        typedef typename"},
{"lineNum":"  107","line":"            spirit::detail::extract_locals<template_params>::type"},
{"lineNum":"  108","line":"        locals_type;"},
{"lineNum":"  109","line":""},
{"lineNum":"  110","line":"        // The rule\'s skip-parser type"},
{"lineNum":"  111","line":"        typedef typename"},
{"lineNum":"  112","line":"            spirit::detail::extract_component<"},
{"lineNum":"  113","line":"                qi::domain, template_params>::type"},
{"lineNum":"  114","line":"        skipper_type;"},
{"lineNum":"  115","line":""},
{"lineNum":"  116","line":"        // The rule\'s encoding type"},
{"lineNum":"  117","line":"        typedef typename"},
{"lineNum":"  118","line":"            spirit::detail::extract_encoding<template_params>::type"},
{"lineNum":"  119","line":"        encoding_type;"},
{"lineNum":"  120","line":""},
{"lineNum":"  121","line":"        // The rule\'s signature"},
{"lineNum":"  122","line":"        typedef typename"},
{"lineNum":"  123","line":"            spirit::detail::extract_sig<template_params, encoding_type, qi::domain>::type"},
{"lineNum":"  124","line":"        sig_type;"},
{"lineNum":"  125","line":""},
{"lineNum":"  126","line":"        // This is the rule\'s attribute type"},
{"lineNum":"  127","line":"        typedef typename"},
{"lineNum":"  128","line":"            spirit::detail::attr_from_sig<sig_type>::type"},
{"lineNum":"  129","line":"        attr_type;"},
{"lineNum":"  130","line":"        typedef typename add_reference<attr_type>::type attr_reference_type;"},
{"lineNum":"  131","line":""},
{"lineNum":"  132","line":"        // parameter_types is a sequence of types passed as parameters to the rule"},
{"lineNum":"  133","line":"        typedef typename"},
{"lineNum":"  134","line":"            spirit::detail::params_from_sig<sig_type>::type"},
{"lineNum":"  135","line":"        parameter_types;"},
{"lineNum":"  136","line":""},
{"lineNum":"  137","line":"        static size_t const params_size ="},
{"lineNum":"  138","line":"            fusion::result_of::size<parameter_types>::type::value;"},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":"        typedef context<"},
{"lineNum":"  141","line":"            fusion::cons<attr_reference_type, parameter_types>"},
{"lineNum":"  142","line":"          , locals_type>"},
{"lineNum":"  143","line":"        context_type;"},
{"lineNum":"  144","line":""},
{"lineNum":"  145","line":"        typedef function<"},
{"lineNum":"  146","line":"            bool(Iterator& first, Iterator const& last"},
{"lineNum":"  147","line":"              , context_type& context"},
{"lineNum":"  148","line":"              , skipper_type const& skipper"},
{"lineNum":"  149","line":"            )>"},
{"lineNum":"  150","line":"        function_type;"},
{"lineNum":"  151","line":""},
{"lineNum":"  152","line":"        typedef typename"},
{"lineNum":"  153","line":"            mpl::if_<"},
{"lineNum":"  154","line":"                is_same<encoding_type, unused_type>"},
{"lineNum":"  155","line":"              , unused_type"},
{"lineNum":"  156","line":"              , tag::char_code<tag::encoding, encoding_type>"},
{"lineNum":"  157","line":"            >::type"},
{"lineNum":"  158","line":"        encoding_modifier_type;"},
{"lineNum":"  159","line":""},
{"lineNum":"  160","line":"        explicit rule(std::string const& name = \"unnamed-rule\")"},
{"lineNum":"  161","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"  162","line":"          , name_(name)","class":"lineNoCov","hits":"0","possible_hits":"28",},
{"lineNum":"  163","line":"        {"},
{"lineNum":"  164","line":"        }"},
{"lineNum":"  165","line":""},
{"lineNum":"  166","line":"        rule(rule const& rhs)"},
{"lineNum":"  167","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"  168","line":"          , name_(rhs.name_)"},
{"lineNum":"  169","line":"          , f(rhs.f)"},
{"lineNum":"  170","line":"        {"},
{"lineNum":"  171","line":"        }"},
{"lineNum":"  172","line":""},
{"lineNum":"  173","line":"        template <typename Auto, typename Expr>"},
{"lineNum":"  174","line":"        static void define(rule& /*lhs*/, Expr const& /*expr*/, mpl::false_)"},
{"lineNum":"  175","line":"        {"},
{"lineNum":"  176","line":"            // Report invalid expression error as early as possible."},
{"lineNum":"  177","line":"            // If you got an error_invalid_expression error message here,"},
{"lineNum":"  178","line":"            // then the expression (expr) is not a valid spirit qi expression."},
{"lineNum":"  179","line":"            BOOST_SPIRIT_ASSERT_MATCH(qi::domain, Expr);"},
{"lineNum":"  180","line":"        }"},
{"lineNum":"  181","line":""},
{"lineNum":"  182","line":"        template <typename Auto, typename Expr>"},
{"lineNum":"  183","line":"        static void define(rule& lhs, Expr const& expr, mpl::true_)"},
{"lineNum":"  184","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"18",},
{"lineNum":"  185","line":"            lhs.f = detail::bind_parser<Auto>(","class":"lineNoCov","hits":"0","possible_hits":"21",},
{"lineNum":"  186","line":"                compile<qi::domain>(expr, encoding_modifier_type()));"},
{"lineNum":"  187","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"18",},
{"lineNum":"  188","line":""},
{"lineNum":"  189","line":"        template <typename Expr>"},
{"lineNum":"  190","line":"        rule(Expr const& expr, std::string const& name = \"unnamed-rule\")","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  191","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"  192","line":"          , name_(name)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  193","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  194","line":"            define<mpl::false_>(*this, expr, traits::matches<qi::domain, Expr>());"},
{"lineNum":"  195","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  196","line":""},
{"lineNum":"  197","line":"        rule& operator=(rule const& rhs)"},
{"lineNum":"  198","line":"        {"},
{"lineNum":"  199","line":"            // The following assertion fires when you try to initialize a rule"},
{"lineNum":"  200","line":"            // from an uninitialized one. Did you mean to refer to the right"},
{"lineNum":"  201","line":"            // hand side rule instead of assigning from it? In this case you"},
{"lineNum":"  202","line":"            // should write lhs = rhs.alias();"},
{"lineNum":"  203","line":"            BOOST_ASSERT(rhs.f && \"Did you mean rhs.alias() instead of rhs?\");"},
{"lineNum":"  204","line":""},
{"lineNum":"  205","line":"            f = rhs.f;"},
{"lineNum":"  206","line":"            name_ = rhs.name_;"},
{"lineNum":"  207","line":"            return *this;"},
{"lineNum":"  208","line":"        }"},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"        std::string const& name() const"},
{"lineNum":"  211","line":"        {"},
{"lineNum":"  212","line":"            return name_;"},
{"lineNum":"  213","line":"        }"},
{"lineNum":"  214","line":""},
{"lineNum":"  215","line":"        void name(std::string const& str)"},
{"lineNum":"  216","line":"        {"},
{"lineNum":"  217","line":"            name_ = str;"},
{"lineNum":"  218","line":"        }"},
{"lineNum":"  219","line":""},
{"lineNum":"  220","line":"        template <typename Expr>"},
{"lineNum":"  221","line":"        rule& operator=(Expr const& expr)"},
{"lineNum":"  222","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  223","line":"            define<mpl::false_>(*this, expr, traits::matches<qi::domain, Expr>());","class":"lineNoCov","hits":"0","possible_hits":"19",},
{"lineNum":"  224","line":"            return *this;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  225","line":"        }"},
{"lineNum":"  226","line":""},
{"lineNum":"  227","line":"// VC7.1 has problems to resolve \'rule\' without explicit template parameters"},
{"lineNum":"  228","line":"#if !BOOST_WORKAROUND(BOOST_MSVC, < 1400)"},
{"lineNum":"  229","line":"        // g++ 3.3 barfs if this is a member function :("},
{"lineNum":"  230","line":"        template <typename Expr>"},
{"lineNum":"  231","line":"        friend rule& operator%=(rule& r, Expr const& expr)"},
{"lineNum":"  232","line":"        {"},
{"lineNum":"  233","line":"            define<mpl::true_>(r, expr, traits::matches<qi::domain, Expr>());"},
{"lineNum":"  234","line":"            return r;"},
{"lineNum":"  235","line":"        }"},
{"lineNum":"  236","line":""},
{"lineNum":"  237","line":"#if defined(BOOST_NO_CXX11_RVALUE_REFERENCES)"},
{"lineNum":"  238","line":"        // non-const version needed to suppress proto\'s %= kicking in"},
{"lineNum":"  239","line":"        template <typename Expr>"},
{"lineNum":"  240","line":"        friend rule& operator%=(rule& r, Expr& expr)"},
{"lineNum":"  241","line":"        {"},
{"lineNum":"  242","line":"            return r %= static_cast<Expr const&>(expr);"},
{"lineNum":"  243","line":"        }"},
{"lineNum":"  244","line":"#else"},
{"lineNum":"  245","line":"        // for rvalue references"},
{"lineNum":"  246","line":"        template <typename Expr>"},
{"lineNum":"  247","line":"        friend rule& operator%=(rule& r, Expr&& expr)"},
{"lineNum":"  248","line":"        {"},
{"lineNum":"  249","line":"            define<mpl::true_>(r, expr, traits::matches<qi::domain, Expr>());"},
{"lineNum":"  250","line":"            return r;"},
{"lineNum":"  251","line":"        }"},
{"lineNum":"  252","line":"#endif"},
{"lineNum":"  253","line":""},
{"lineNum":"  254","line":"#else"},
{"lineNum":"  255","line":"        // both friend functions have to be defined out of class as VC7.1"},
{"lineNum":"  256","line":"        // will complain otherwise"},
{"lineNum":"  257","line":"        template <typename OutputIterator_, typename T1_, typename T2_"},
{"lineNum":"  258","line":"          , typename T3_, typename T4_, typename Expr>"},
{"lineNum":"  259","line":"        friend rule<OutputIterator_, T1_, T2_, T3_, T4_>& operator%=("},
{"lineNum":"  260","line":"            rule<OutputIterator_, T1_, T2_, T3_, T4_>& r, Expr const& expr);"},
{"lineNum":"  261","line":""},
{"lineNum":"  262","line":"        // non-const version needed to suppress proto\'s %= kicking in"},
{"lineNum":"  263","line":"        template <typename OutputIterator_, typename T1_, typename T2_"},
{"lineNum":"  264","line":"          , typename T3_, typename T4_, typename Expr>"},
{"lineNum":"  265","line":"        friend rule<OutputIterator_, T1_, T2_, T3_, T4_>& operator%=("},
{"lineNum":"  266","line":"            rule<OutputIterator_, T1_, T2_, T3_, T4_>& r, Expr& expr);"},
{"lineNum":"  267","line":"#endif"},
{"lineNum":"  268","line":""},
{"lineNum":"  269","line":"        template <typename Context, typename Iterator_>"},
{"lineNum":"  270","line":"        struct attribute"},
{"lineNum":"  271","line":"        {"},
{"lineNum":"  272","line":"            typedef attr_type type;"},
{"lineNum":"  273","line":"        };"},
{"lineNum":"  274","line":""},
{"lineNum":"  275","line":"        template <typename Context, typename Skipper, typename Attribute>"},
{"lineNum":"  276","line":"        bool parse(Iterator& first, Iterator const& last"},
{"lineNum":"  277","line":"          , Context& /*context*/, Skipper const& skipper"},
{"lineNum":"  278","line":"          , Attribute& attr_param) const"},
{"lineNum":"  279","line":"        {"},
{"lineNum":"  280","line":"            BOOST_STATIC_ASSERT_MSG((is_same<skipper_type, unused_type>::value ||"},
{"lineNum":"  281","line":"                !is_same<Skipper, unused_type>::value),"},
{"lineNum":"  282","line":"                \"The rule was instantiated with a skipper type but you have not pass any. \""},
{"lineNum":"  283","line":"                \"Did you use `parse` instead of `phrase_parse`?\");"},
{"lineNum":"  284","line":"            BOOST_STATIC_ASSERT_MSG("},
{"lineNum":"  285","line":"                (is_convertible<Skipper const&, skipper_type>::value),"},
{"lineNum":"  286","line":"                \"The passed skipper is not compatible/convertible to one \""},
{"lineNum":"  287","line":"                \"that the rule was instantiated with\");"},
{"lineNum":"  288","line":"            if (f)","class":"lineNoCov","hits":"0","possible_hits":"93",},
{"lineNum":"  289","line":"            {"},
{"lineNum":"  290","line":"                // do a preskip if this is an implied lexeme"},
{"lineNum":"  291","line":"                if (is_same<skipper_type, unused_type>::value)"},
{"lineNum":"  292","line":"                    qi::skip_over(first, last, skipper);"},
{"lineNum":"  293","line":""},
{"lineNum":"  294","line":"                // do down-stream transformation, provides attribute for"},
{"lineNum":"  295","line":"                // rhs parser"},
{"lineNum":"  296","line":"                typedef traits::transform_attribute<"},
{"lineNum":"  297","line":"                    Attribute, attr_type, domain>"},
{"lineNum":"  298","line":"                transform;"},
{"lineNum":"  299","line":""},
{"lineNum":"  300","line":"                typename transform::type attr_ = transform::pre(attr_param);"},
{"lineNum":"  301","line":""},
{"lineNum":"  302","line":"                // If you are seeing a compilation error here, you are probably"},
{"lineNum":"  303","line":"                // trying to use a rule or a grammar which has inherited"},
{"lineNum":"  304","line":"                // attributes, without passing values for them."},
{"lineNum":"  305","line":"                context_type context(attr_);"},
{"lineNum":"  306","line":""},
{"lineNum":"  307","line":"                // If you are seeing a compilation error here stating that the"},
{"lineNum":"  308","line":"                // fourth parameter can\'t be converted to a required target type"},
{"lineNum":"  309","line":"                // then you are probably trying to use a rule or a grammar with"},
{"lineNum":"  310","line":"                // an incompatible skipper type."},
{"lineNum":"  311","line":"                if (f(first, last, context, skipper))"},
{"lineNum":"  312","line":"                {"},
{"lineNum":"  313","line":"                    // do up-stream transformation, this integrates the results"},
{"lineNum":"  314","line":"                    // back into the original attribute value, if appropriate"},
{"lineNum":"  315","line":"                    transform::post(attr_param, attr_);"},
{"lineNum":"  316","line":"                    return true;"},
{"lineNum":"  317","line":"                }"},
{"lineNum":"  318","line":""},
{"lineNum":"  319","line":"                // inform attribute transformation of failed rhs"},
{"lineNum":"  320","line":"                transform::fail(attr_param);"},
{"lineNum":"  321","line":"            }"},
{"lineNum":"  322","line":"            return false;"},
{"lineNum":"  323","line":"        }"},
{"lineNum":"  324","line":""},
{"lineNum":"  325","line":"        template <typename Context, typename Skipper"},
{"lineNum":"  326","line":"          , typename Attribute, typename Params>"},
{"lineNum":"  327","line":"        bool parse(Iterator& first, Iterator const& last"},
{"lineNum":"  328","line":"          , Context& caller_context, Skipper const& skipper"},
{"lineNum":"  329","line":"          , Attribute& attr_param, Params const& params) const"},
{"lineNum":"  330","line":"        {"},
{"lineNum":"  331","line":"            BOOST_STATIC_ASSERT_MSG((is_same<skipper_type, unused_type>::value ||"},
{"lineNum":"  332","line":"                !is_same<Skipper, unused_type>::value),"},
{"lineNum":"  333","line":"                \"The rule was instantiated with a skipper type but you have not pass any. \""},
{"lineNum":"  334","line":"                \"Did you use `parse` instead of `phrase_parse`?\");"},
{"lineNum":"  335","line":"            BOOST_STATIC_ASSERT_MSG("},
{"lineNum":"  336","line":"                (is_convertible<Skipper const&, skipper_type>::value),"},
{"lineNum":"  337","line":"                \"The passed skipper is not compatible/convertible to one \""},
{"lineNum":"  338","line":"                \"that the rule was instantiated with\");"},
{"lineNum":"  339","line":"            if (f)"},
{"lineNum":"  340","line":"            {"},
{"lineNum":"  341","line":"                // do a preskip if this is an implied lexeme"},
{"lineNum":"  342","line":"                if (is_same<skipper_type, unused_type>::value)"},
{"lineNum":"  343","line":"                    qi::skip_over(first, last, skipper);"},
{"lineNum":"  344","line":""},
{"lineNum":"  345","line":"                // do down-stream transformation, provides attribute for"},
{"lineNum":"  346","line":"                // rhs parser"},
{"lineNum":"  347","line":"                typedef traits::transform_attribute<"},
{"lineNum":"  348","line":"                    Attribute, attr_type, domain>"},
{"lineNum":"  349","line":"                transform;"},
{"lineNum":"  350","line":""},
{"lineNum":"  351","line":"                typename transform::type attr_ = transform::pre(attr_param);"},
{"lineNum":"  352","line":""},
{"lineNum":"  353","line":"                // If you are seeing a compilation error here, you are probably"},
{"lineNum":"  354","line":"                // trying to use a rule or a grammar which has inherited"},
{"lineNum":"  355","line":"                // attributes, passing values of incompatible types for them."},
{"lineNum":"  356","line":"                context_type context(attr_, params, caller_context);"},
{"lineNum":"  357","line":""},
{"lineNum":"  358","line":"                // If you are seeing a compilation error here stating that the"},
{"lineNum":"  359","line":"                // fourth parameter can\'t be converted to a required target type"},
{"lineNum":"  360","line":"                // then you are probably trying to use a rule or a grammar with"},
{"lineNum":"  361","line":"                // an incompatible skipper type."},
{"lineNum":"  362","line":"                if (f(first, last, context, skipper))"},
{"lineNum":"  363","line":"                {"},
{"lineNum":"  364","line":"                    // do up-stream transformation, this integrates the results"},
{"lineNum":"  365","line":"                    // back into the original attribute value, if appropriate"},
{"lineNum":"  366","line":"                    transform::post(attr_param, attr_);"},
{"lineNum":"  367","line":"                    return true;"},
{"lineNum":"  368","line":"                }"},
{"lineNum":"  369","line":""},
{"lineNum":"  370","line":"                // inform attribute transformation of failed rhs"},
{"lineNum":"  371","line":"                transform::fail(attr_param);"},
{"lineNum":"  372","line":"            }"},
{"lineNum":"  373","line":"            return false;"},
{"lineNum":"  374","line":"        }"},
{"lineNum":"  375","line":""},
{"lineNum":"  376","line":"        template <typename Context>"},
{"lineNum":"  377","line":"        info what(Context& /*context*/) const"},
{"lineNum":"  378","line":"        {"},
{"lineNum":"  379","line":"            return info(name_);"},
{"lineNum":"  380","line":"        }"},
{"lineNum":"  381","line":""},
{"lineNum":"  382","line":"        reference_ alias() const"},
{"lineNum":"  383","line":"        {"},
{"lineNum":"  384","line":"            return reference_(*this);"},
{"lineNum":"  385","line":"        }"},
{"lineNum":"  386","line":""},
{"lineNum":"  387","line":"        typename proto::terminal<this_type>::type copy() const"},
{"lineNum":"  388","line":"        {"},
{"lineNum":"  389","line":"            typename proto::terminal<this_type>::type result = {*this};"},
{"lineNum":"  390","line":"            return result;"},
{"lineNum":"  391","line":"        }"},
{"lineNum":"  392","line":""},
{"lineNum":"  393","line":"        // bring in the operator() overloads"},
{"lineNum":"  394","line":"        rule const& get_parameterized_subject() const { return *this; }"},
{"lineNum":"  395","line":"        typedef rule parameterized_subject_type;"},
{"lineNum":"  396","line":"        #include <boost/spirit/home/qi/nonterminal/detail/fcall.hpp>"},
{"lineNum":"  397","line":""},
{"lineNum":"  398","line":"        std::string name_;"},
{"lineNum":"  399","line":"        function_type f;"},
{"lineNum":"  400","line":"    };"},
{"lineNum":"  401","line":""},
{"lineNum":"  402","line":"#if BOOST_WORKAROUND(BOOST_MSVC, < 1400)"},
{"lineNum":"  403","line":"    template <typename OutputIterator_, typename T1_, typename T2_"},
{"lineNum":"  404","line":"      , typename T3_, typename T4_, typename Expr>"},
{"lineNum":"  405","line":"    rule<OutputIterator_, T1_, T2_, T3_, T4_>& operator%=("},
{"lineNum":"  406","line":"        rule<OutputIterator_, T1_, T2_, T3_, T4_>& r, Expr const& expr)"},
{"lineNum":"  407","line":"    {"},
{"lineNum":"  408","line":"        // Report invalid expression error as early as possible."},
{"lineNum":"  409","line":"        // If you got an error_invalid_expression error message here,"},
{"lineNum":"  410","line":"        // then the expression (expr) is not a valid spirit qi expression."},
{"lineNum":"  411","line":"        BOOST_SPIRIT_ASSERT_MATCH(qi::domain, Expr);"},
{"lineNum":"  412","line":""},
{"lineNum":"  413","line":"        typedef typename"},
{"lineNum":"  414","line":"            rule<OutputIterator_, T1_, T2_, T3_, T4_>::encoding_modifier_type"},
{"lineNum":"  415","line":"        encoding_modifier_type;"},
{"lineNum":"  416","line":""},
{"lineNum":"  417","line":"        r.f = detail::bind_parser<mpl::true_>("},
{"lineNum":"  418","line":"            compile<qi::domain>(expr, encoding_modifier_type()));"},
{"lineNum":"  419","line":"        return r;"},
{"lineNum":"  420","line":"    }"},
{"lineNum":"  421","line":""},
{"lineNum":"  422","line":"    template <typename Iterator_, typename T1_, typename T2_"},
{"lineNum":"  423","line":"      , typename T3_, typename T4_, typename Expr>"},
{"lineNum":"  424","line":"    rule<Iterator_, T1_, T2_, T3_, T4_>& operator%=("},
{"lineNum":"  425","line":"        rule<Iterator_, T1_, T2_, T3_, T4_>& r, Expr& expr)"},
{"lineNum":"  426","line":"    {"},
{"lineNum":"  427","line":"        return r %= static_cast<Expr const&>(expr);"},
{"lineNum":"  428","line":"    }"},
{"lineNum":"  429","line":"#endif"},
{"lineNum":"  430","line":"}}}"},
{"lineNum":"  431","line":""},
{"lineNum":"  432","line":"namespace boost { namespace spirit { namespace traits"},
{"lineNum":"  433","line":"{"},
{"lineNum":"  434","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  435","line":"    template <"},
{"lineNum":"  436","line":"        typename IteratorA, typename IteratorB, typename Attribute"},
{"lineNum":"  437","line":"      , typename Context, typename T1, typename T2, typename T3, typename T4>"},
{"lineNum":"  438","line":"    struct handles_container<"},
{"lineNum":"  439","line":"        qi::rule<IteratorA, T1, T2, T3, T4>, Attribute, Context, IteratorB>"},
{"lineNum":"  440","line":"      : traits::is_container<"},
{"lineNum":"  441","line":"          typename attribute_of<"},
{"lineNum":"  442","line":"              qi::rule<IteratorA, T1, T2, T3, T4>, Context, IteratorB"},
{"lineNum":"  443","line":"          >::type"},
{"lineNum":"  444","line":"        >"},
{"lineNum":"  445","line":"    {};"},
{"lineNum":"  446","line":"}}}"},
{"lineNum":"  447","line":""},
{"lineNum":"  448","line":"#if defined(BOOST_MSVC)"},
{"lineNum":"  449","line":"# pragma warning(pop)"},
{"lineNum":"  450","line":"#endif"},
{"lineNum":"  451","line":""},
{"lineNum":"  452","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 13, "covered" : 0,};
var merged_data = [];
