var data = {lines:[
{"lineNum":"    1","line":"/*============================================================================="},
{"lineNum":"    2","line":"    Copyright (c) 2001-2011 Joel de Guzman"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":"    Distributed under the Boost Software License, Version 1.0. (See accompanying"},
{"lineNum":"    5","line":"    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    6","line":"==============================================================================*/"},
{"lineNum":"    7","line":"#if !defined(BOOST_SPIRIT_SYMBOLS_MARCH_11_2007_1055AM)"},
{"lineNum":"    8","line":"#define BOOST_SPIRIT_SYMBOLS_MARCH_11_2007_1055AM"},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#if defined(_MSC_VER)"},
{"lineNum":"   11","line":"#pragma once"},
{"lineNum":"   12","line":"#endif"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"#include <boost/spirit/home/qi/domain.hpp>"},
{"lineNum":"   15","line":"#include <boost/spirit/home/qi/skip_over.hpp>"},
{"lineNum":"   16","line":"#include <boost/spirit/home/qi/string/tst.hpp>"},
{"lineNum":"   17","line":"#include <boost/spirit/home/qi/reference.hpp>"},
{"lineNum":"   18","line":"#include <boost/spirit/home/qi/meta_compiler.hpp>"},
{"lineNum":"   19","line":"#include <boost/spirit/home/qi/detail/assign_to.hpp>"},
{"lineNum":"   20","line":"#include <boost/spirit/home/qi/parser.hpp>"},
{"lineNum":"   21","line":"#include <boost/spirit/home/support/detail/get_encoding.hpp>"},
{"lineNum":"   22","line":"#include <boost/spirit/home/support/modify.hpp>"},
{"lineNum":"   23","line":"#include <boost/spirit/home/support/info.hpp>"},
{"lineNum":"   24","line":"#include <boost/spirit/home/support/unused.hpp>"},
{"lineNum":"   25","line":"#include <boost/spirit/home/support/string_traits.hpp>"},
{"lineNum":"   26","line":""},
{"lineNum":"   27","line":"#include <boost/range/begin.hpp>"},
{"lineNum":"   28","line":"#include <boost/range/end.hpp>"},
{"lineNum":"   29","line":"#include <boost/shared_ptr.hpp>"},
{"lineNum":"   30","line":""},
{"lineNum":"   31","line":"#if defined(BOOST_MSVC)"},
{"lineNum":"   32","line":"# pragma warning(push)"},
{"lineNum":"   33","line":"# pragma warning(disable: 4355) // \'this\' : used in base member initializer list warning"},
{"lineNum":"   34","line":"#endif"},
{"lineNum":"   35","line":""},
{"lineNum":"   36","line":"namespace boost { namespace spirit { namespace qi"},
{"lineNum":"   37","line":"{"},
{"lineNum":"   38","line":"    template <"},
{"lineNum":"   39","line":"        typename Char = char"},
{"lineNum":"   40","line":"      , typename T = unused_type"},
{"lineNum":"   41","line":"      , typename Lookup = tst<Char, T>"},
{"lineNum":"   42","line":"      , typename Filter = tst_pass_through>"},
{"lineNum":"   43","line":"    struct symbols","class":"lineNoCov","hits":"0","possible_hits":"27",},
{"lineNum":"   44","line":"      : proto::extends<"},
{"lineNum":"   45","line":"            typename proto::terminal<"},
{"lineNum":"   46","line":"                reference<symbols<Char, T, Lookup, Filter> >"},
{"lineNum":"   47","line":"            >::type"},
{"lineNum":"   48","line":"          , symbols<Char, T, Lookup, Filter>"},
{"lineNum":"   49","line":"        >"},
{"lineNum":"   50","line":"      , primitive_parser<symbols<Char, T, Lookup, Filter> >"},
{"lineNum":"   51","line":"    {"},
{"lineNum":"   52","line":"        typedef Char char_type; // the character type"},
{"lineNum":"   53","line":"        typedef T value_type; // the value associated with each entry"},
{"lineNum":"   54","line":"        typedef symbols<Char, T, Lookup, Filter> this_type;"},
{"lineNum":"   55","line":"        typedef reference<this_type> reference_;"},
{"lineNum":"   56","line":"        typedef typename proto::terminal<reference_>::type terminal;"},
{"lineNum":"   57","line":"        typedef proto::extends<terminal, this_type> base_type;"},
{"lineNum":"   58","line":""},
{"lineNum":"   59","line":"        template <typename Context, typename Iterator>"},
{"lineNum":"   60","line":"        struct attribute"},
{"lineNum":"   61","line":"        {"},
{"lineNum":"   62","line":"            typedef value_type type;"},
{"lineNum":"   63","line":"        };"},
{"lineNum":"   64","line":""},
{"lineNum":"   65","line":"        symbols(std::string const& name = \"symbols\")"},
{"lineNum":"   66","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"   67","line":"          , add(*this)"},
{"lineNum":"   68","line":"          , remove(*this)"},
{"lineNum":"   69","line":"          , lookup(new Lookup())","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   70","line":"          , name_(name)","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   71","line":"        {"},
{"lineNum":"   72","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"11",},
{"lineNum":"   73","line":""},
{"lineNum":"   74","line":"        symbols(symbols const& syms)"},
{"lineNum":"   75","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"   76","line":"          , add(*this)"},
{"lineNum":"   77","line":"          , remove(*this)"},
{"lineNum":"   78","line":"          , lookup(syms.lookup)"},
{"lineNum":"   79","line":"          , name_(syms.name_)","class":"lineNoCov","hits":"0","possible_hits":"48",},
{"lineNum":"   80","line":"        {"},
{"lineNum":"   81","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"        template <typename Filter_>"},
{"lineNum":"   84","line":"        symbols(symbols<Char, T, Lookup, Filter_> const& syms)"},
{"lineNum":"   85","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"   86","line":"          , add(*this)"},
{"lineNum":"   87","line":"          , remove(*this)"},
{"lineNum":"   88","line":"          , lookup(syms.lookup)"},
{"lineNum":"   89","line":"          , name_(syms.name_)","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"   90","line":"        {"},
{"lineNum":"   91","line":"        }"},
{"lineNum":"   92","line":""},
{"lineNum":"   93","line":"        template <typename Symbols>"},
{"lineNum":"   94","line":"        symbols(Symbols const& syms, std::string const& name = \"symbols\")"},
{"lineNum":"   95","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"   96","line":"          , add(*this)"},
{"lineNum":"   97","line":"          , remove(*this)"},
{"lineNum":"   98","line":"          , lookup(new Lookup())"},
{"lineNum":"   99","line":"          , name_(name)"},
{"lineNum":"  100","line":"        {"},
{"lineNum":"  101","line":"            typename range_const_iterator<Symbols>::type si = boost::begin(syms);"},
{"lineNum":"  102","line":"            while (si != boost::end(syms))"},
{"lineNum":"  103","line":"                add(*si++);"},
{"lineNum":"  104","line":"        }"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"        template <typename Symbols, typename Data>"},
{"lineNum":"  107","line":"        symbols(Symbols const& syms, Data const& data"},
{"lineNum":"  108","line":"              , std::string const& name = \"symbols\")"},
{"lineNum":"  109","line":"          : base_type(terminal::make(reference_(*this)))"},
{"lineNum":"  110","line":"          , add(*this)"},
{"lineNum":"  111","line":"          , remove(*this)"},
{"lineNum":"  112","line":"          , lookup(new Lookup())"},
{"lineNum":"  113","line":"          , name_(name)"},
{"lineNum":"  114","line":"        {"},
{"lineNum":"  115","line":"            typename range_const_iterator<Symbols>::type si = boost::begin(syms);"},
{"lineNum":"  116","line":"            typename range_const_iterator<Data>::type di = boost::begin(data);"},
{"lineNum":"  117","line":"            while (si != boost::end(syms))"},
{"lineNum":"  118","line":"                add(*si++, *di++);"},
{"lineNum":"  119","line":"        }"},
{"lineNum":"  120","line":""},
{"lineNum":"  121","line":"        symbols&"},
{"lineNum":"  122","line":"        operator=(symbols const& rhs)"},
{"lineNum":"  123","line":"        {"},
{"lineNum":"  124","line":"            name_ = rhs.name_;"},
{"lineNum":"  125","line":"            *lookup = *rhs.lookup;"},
{"lineNum":"  126","line":"            return *this;"},
{"lineNum":"  127","line":"        }"},
{"lineNum":"  128","line":""},
{"lineNum":"  129","line":"        template <typename Filter_>"},
{"lineNum":"  130","line":"        symbols&"},
{"lineNum":"  131","line":"        operator=(symbols<Char, T, Lookup, Filter_> const& rhs)"},
{"lineNum":"  132","line":"        {"},
{"lineNum":"  133","line":"            name_ = rhs.name_;"},
{"lineNum":"  134","line":"            *lookup = *rhs.lookup;"},
{"lineNum":"  135","line":"            return *this;"},
{"lineNum":"  136","line":"        }"},
{"lineNum":"  137","line":""},
{"lineNum":"  138","line":"        void clear()"},
{"lineNum":"  139","line":"        {"},
{"lineNum":"  140","line":"            lookup->clear();"},
{"lineNum":"  141","line":"        }"},
{"lineNum":"  142","line":""},
{"lineNum":"  143","line":"        struct adder;"},
{"lineNum":"  144","line":"        struct remover;"},
{"lineNum":"  145","line":""},
{"lineNum":"  146","line":"        template <typename Str>"},
{"lineNum":"  147","line":"        adder const&"},
{"lineNum":"  148","line":"        operator=(Str const& str)"},
{"lineNum":"  149","line":"        {"},
{"lineNum":"  150","line":"            lookup->clear();"},
{"lineNum":"  151","line":"            return add(str);"},
{"lineNum":"  152","line":"        }"},
{"lineNum":"  153","line":""},
{"lineNum":"  154","line":"        template <typename Str>"},
{"lineNum":"  155","line":"        friend adder const&"},
{"lineNum":"  156","line":"        operator+=(symbols& sym, Str const& str)"},
{"lineNum":"  157","line":"        {"},
{"lineNum":"  158","line":"            return sym.add(str);"},
{"lineNum":"  159","line":"        }"},
{"lineNum":"  160","line":""},
{"lineNum":"  161","line":"        template <typename Str>"},
{"lineNum":"  162","line":"        friend remover const&"},
{"lineNum":"  163","line":"        operator-=(symbols& sym, Str const& str)"},
{"lineNum":"  164","line":"        {"},
{"lineNum":"  165","line":"            return sym.remove(str);"},
{"lineNum":"  166","line":"        }"},
{"lineNum":"  167","line":""},
{"lineNum":"  168","line":"#if defined(BOOST_NO_CXX11_RVALUE_REFERENCES)"},
{"lineNum":"  169","line":"        // non-const version needed to suppress proto\'s += kicking in"},
{"lineNum":"  170","line":"        template <typename Str>"},
{"lineNum":"  171","line":"        friend adder const&"},
{"lineNum":"  172","line":"        operator+=(symbols& sym, Str& str)"},
{"lineNum":"  173","line":"        {"},
{"lineNum":"  174","line":"            return sym.add(str);"},
{"lineNum":"  175","line":"        }"},
{"lineNum":"  176","line":""},
{"lineNum":"  177","line":"        // non-const version needed to suppress proto\'s -= kicking in"},
{"lineNum":"  178","line":"        template <typename Str>"},
{"lineNum":"  179","line":"        friend remover const&"},
{"lineNum":"  180","line":"        operator-=(symbols& sym, Str& str)"},
{"lineNum":"  181","line":"        {"},
{"lineNum":"  182","line":"            return sym.remove(str);"},
{"lineNum":"  183","line":"        }"},
{"lineNum":"  184","line":"#else"},
{"lineNum":"  185","line":"        // for rvalue references"},
{"lineNum":"  186","line":"        template <typename Str>"},
{"lineNum":"  187","line":"        friend adder const&"},
{"lineNum":"  188","line":"        operator+=(symbols& sym, Str&& str)"},
{"lineNum":"  189","line":"        {"},
{"lineNum":"  190","line":"            return sym.add(str);"},
{"lineNum":"  191","line":"        }"},
{"lineNum":"  192","line":""},
{"lineNum":"  193","line":"        // for rvalue references"},
{"lineNum":"  194","line":"        template <typename Str>"},
{"lineNum":"  195","line":"        friend remover const&"},
{"lineNum":"  196","line":"        operator-=(symbols& sym, Str&& str)"},
{"lineNum":"  197","line":"        {"},
{"lineNum":"  198","line":"            return sym.remove(str);"},
{"lineNum":"  199","line":"        }"},
{"lineNum":"  200","line":"#endif"},
{"lineNum":"  201","line":"        template <typename F>"},
{"lineNum":"  202","line":"        void for_each(F f) const"},
{"lineNum":"  203","line":"        {"},
{"lineNum":"  204","line":"            lookup->for_each(f);"},
{"lineNum":"  205","line":"        }"},
{"lineNum":"  206","line":""},
{"lineNum":"  207","line":"        template <typename Str>"},
{"lineNum":"  208","line":"        value_type& at(Str const& str)"},
{"lineNum":"  209","line":"        {"},
{"lineNum":"  210","line":"            return *lookup->add(traits::get_begin<Char>(str)"},
{"lineNum":"  211","line":"                , traits::get_end<Char>(str), T());"},
{"lineNum":"  212","line":"        }"},
{"lineNum":"  213","line":""},
{"lineNum":"  214","line":"        template <typename Iterator>"},
{"lineNum":"  215","line":"        value_type* prefix_find(Iterator& first, Iterator const& last)"},
{"lineNum":"  216","line":"        {"},
{"lineNum":"  217","line":"            return lookup->find(first, last, Filter());"},
{"lineNum":"  218","line":"        }"},
{"lineNum":"  219","line":""},
{"lineNum":"  220","line":"        template <typename Iterator>"},
{"lineNum":"  221","line":"        value_type const* prefix_find(Iterator& first, Iterator const& last) const"},
{"lineNum":"  222","line":"        {"},
{"lineNum":"  223","line":"            return lookup->find(first, last, Filter());"},
{"lineNum":"  224","line":"        }"},
{"lineNum":"  225","line":""},
{"lineNum":"  226","line":"        template <typename Str>"},
{"lineNum":"  227","line":"        value_type* find(Str const& str)"},
{"lineNum":"  228","line":"        {"},
{"lineNum":"  229","line":"            return find_impl(traits::get_begin<Char>(str)"},
{"lineNum":"  230","line":"                , traits::get_end<Char>(str));"},
{"lineNum":"  231","line":"        }"},
{"lineNum":"  232","line":""},
{"lineNum":"  233","line":"        template <typename Str>"},
{"lineNum":"  234","line":"        value_type const* find(Str const& str) const"},
{"lineNum":"  235","line":"        {"},
{"lineNum":"  236","line":"            return find_impl(traits::get_begin<Char>(str)"},
{"lineNum":"  237","line":"                , traits::get_end<Char>(str));"},
{"lineNum":"  238","line":"        }"},
{"lineNum":"  239","line":""},
{"lineNum":"  240","line":"private:"},
{"lineNum":"  241","line":"        template <typename Iterator>"},
{"lineNum":"  242","line":"        value_type* find_impl(Iterator begin, Iterator end)"},
{"lineNum":"  243","line":"        {"},
{"lineNum":"  244","line":"            value_type* r = lookup->find(begin, end, Filter());"},
{"lineNum":"  245","line":"            return begin == end ? r : 0;"},
{"lineNum":"  246","line":"        }"},
{"lineNum":"  247","line":""},
{"lineNum":"  248","line":"        template <typename Iterator>"},
{"lineNum":"  249","line":"        value_type const* find_impl(Iterator begin, Iterator end) const"},
{"lineNum":"  250","line":"        {"},
{"lineNum":"  251","line":"            value_type const* r = lookup->find(begin, end, Filter());"},
{"lineNum":"  252","line":"            return begin == end ? r : 0;"},
{"lineNum":"  253","line":"        }"},
{"lineNum":"  254","line":""},
{"lineNum":"  255","line":"public:"},
{"lineNum":"  256","line":"        template <typename Iterator, typename Context"},
{"lineNum":"  257","line":"          , typename Skipper, typename Attribute>"},
{"lineNum":"  258","line":"        bool parse(Iterator& first, Iterator const& last"},
{"lineNum":"  259","line":"          , Context& /*context*/, Skipper const& skipper, Attribute& attr_) const"},
{"lineNum":"  260","line":"        {"},
{"lineNum":"  261","line":"            qi::skip_over(first, last, skipper);"},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"            if (value_type* val_ptr"},
{"lineNum":"  264","line":"                = lookup->find(first, last, Filter()))","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"  265","line":"            {"},
{"lineNum":"  266","line":"                spirit::traits::assign_to(*val_ptr, attr_);"},
{"lineNum":"  267","line":"                return true;"},
{"lineNum":"  268","line":"            }"},
{"lineNum":"  269","line":"            return false;"},
{"lineNum":"  270","line":"        }"},
{"lineNum":"  271","line":""},
{"lineNum":"  272","line":"        template <typename Context>"},
{"lineNum":"  273","line":"        info what(Context& /*context*/) const"},
{"lineNum":"  274","line":"        {"},
{"lineNum":"  275","line":"            return info(name_);"},
{"lineNum":"  276","line":"        }"},
{"lineNum":"  277","line":""},
{"lineNum":"  278","line":"        void name(std::string const &str)"},
{"lineNum":"  279","line":"        {"},
{"lineNum":"  280","line":"            name_ = str;"},
{"lineNum":"  281","line":"        }"},
{"lineNum":"  282","line":"        std::string const &name() const"},
{"lineNum":"  283","line":"        {"},
{"lineNum":"  284","line":"            return name_;"},
{"lineNum":"  285","line":"        }"},
{"lineNum":"  286","line":""},
{"lineNum":"  287","line":"        struct adder"},
{"lineNum":"  288","line":"        {"},
{"lineNum":"  289","line":"            template <typename, typename = unused_type, typename = unused_type>"},
{"lineNum":"  290","line":"            struct result { typedef adder const& type; };"},
{"lineNum":"  291","line":""},
{"lineNum":"  292","line":"            adder(symbols& sym_)"},
{"lineNum":"  293","line":"              : sym(sym_)","class":"lineNoCov","hits":"0","possible_hits":"66",},
{"lineNum":"  294","line":"            {"},
{"lineNum":"  295","line":"            }"},
{"lineNum":"  296","line":""},
{"lineNum":"  297","line":"            template <typename Iterator>"},
{"lineNum":"  298","line":"            adder const&"},
{"lineNum":"  299","line":"            operator()(Iterator const& first, Iterator const& last, T const& val) const"},
{"lineNum":"  300","line":"            {"},
{"lineNum":"  301","line":"                sym.lookup->add(first, last, val);"},
{"lineNum":"  302","line":"                return *this;"},
{"lineNum":"  303","line":"            }"},
{"lineNum":"  304","line":""},
{"lineNum":"  305","line":"            template <typename Str>"},
{"lineNum":"  306","line":"            adder const&"},
{"lineNum":"  307","line":"            operator()(Str const& s, T const& val = T()) const"},
{"lineNum":"  308","line":"            {","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  309","line":"                sym.lookup->add(traits::get_begin<Char>(s)","class":"lineNoCov","hits":"0","possible_hits":"53",},
{"lineNum":"  310","line":"                  , traits::get_end<Char>(s), val);"},
{"lineNum":"  311","line":"                return *this;","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  312","line":"            }"},
{"lineNum":"  313","line":""},
{"lineNum":"  314","line":"            template <typename Str>"},
{"lineNum":"  315","line":"            adder const&"},
{"lineNum":"  316","line":"            operator,(Str const& s) const"},
{"lineNum":"  317","line":"            {"},
{"lineNum":"  318","line":"                sym.lookup->add(traits::get_begin<Char>(s)"},
{"lineNum":"  319","line":"                  , traits::get_end<Char>(s), T());"},
{"lineNum":"  320","line":"                return *this;"},
{"lineNum":"  321","line":"            }"},
{"lineNum":"  322","line":""},
{"lineNum":"  323","line":"            symbols& sym;"},
{"lineNum":"  324","line":""},
{"lineNum":"  325","line":"        private:"},
{"lineNum":"  326","line":"            // silence MSVC warning C4512: assignment operator could not be generated"},
{"lineNum":"  327","line":"            adder& operator= (adder const&);"},
{"lineNum":"  328","line":"        };"},
{"lineNum":"  329","line":""},
{"lineNum":"  330","line":"        struct remover"},
{"lineNum":"  331","line":"        {"},
{"lineNum":"  332","line":"            template <typename, typename = unused_type, typename = unused_type>"},
{"lineNum":"  333","line":"            struct result { typedef remover const& type; };"},
{"lineNum":"  334","line":""},
{"lineNum":"  335","line":"            remover(symbols& sym_)"},
{"lineNum":"  336","line":"              : sym(sym_)","class":"lineNoCov","hits":"0","possible_hits":"66",},
{"lineNum":"  337","line":"            {"},
{"lineNum":"  338","line":"            }"},
{"lineNum":"  339","line":""},
{"lineNum":"  340","line":"            template <typename Iterator>"},
{"lineNum":"  341","line":"            remover const&"},
{"lineNum":"  342","line":"            operator()(Iterator const& first, Iterator const& last) const"},
{"lineNum":"  343","line":"            {"},
{"lineNum":"  344","line":"                sym.lookup->remove(first, last);"},
{"lineNum":"  345","line":"                return *this;"},
{"lineNum":"  346","line":"            }"},
{"lineNum":"  347","line":""},
{"lineNum":"  348","line":"            template <typename Str>"},
{"lineNum":"  349","line":"            remover const&"},
{"lineNum":"  350","line":"            operator()(Str const& s) const"},
{"lineNum":"  351","line":"            {"},
{"lineNum":"  352","line":"                sym.lookup->remove(traits::get_begin<Char>(s)"},
{"lineNum":"  353","line":"                  , traits::get_end<Char>(s));"},
{"lineNum":"  354","line":"                return *this;"},
{"lineNum":"  355","line":"            }"},
{"lineNum":"  356","line":""},
{"lineNum":"  357","line":"            template <typename Str>"},
{"lineNum":"  358","line":"            remover const&"},
{"lineNum":"  359","line":"            operator,(Str const& s) const"},
{"lineNum":"  360","line":"            {"},
{"lineNum":"  361","line":"                sym.lookup->remove(traits::get_begin<Char>(s)"},
{"lineNum":"  362","line":"                  , traits::get_end<Char>(s));"},
{"lineNum":"  363","line":"                return *this;"},
{"lineNum":"  364","line":"            }"},
{"lineNum":"  365","line":""},
{"lineNum":"  366","line":"            symbols& sym;"},
{"lineNum":"  367","line":""},
{"lineNum":"  368","line":"        private:"},
{"lineNum":"  369","line":"            // silence MSVC warning C4512: assignment operator could not be generated"},
{"lineNum":"  370","line":"            remover& operator= (remover const&);"},
{"lineNum":"  371","line":"        };"},
{"lineNum":"  372","line":""},
{"lineNum":"  373","line":"        adder add;"},
{"lineNum":"  374","line":"        remover remove;"},
{"lineNum":"  375","line":"        shared_ptr<Lookup> lookup;"},
{"lineNum":"  376","line":"        std::string name_;"},
{"lineNum":"  377","line":"    };"},
{"lineNum":"  378","line":""},
{"lineNum":"  379","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  380","line":"    // Parser generators: make_xxx function (objects)"},
{"lineNum":"  381","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  382","line":"    template <typename Char, typename T, typename Lookup"},
{"lineNum":"  383","line":"      , typename Filter, typename Modifiers>"},
{"lineNum":"  384","line":"    struct make_primitive<reference<symbols<Char, T, Lookup, Filter> >, Modifiers>"},
{"lineNum":"  385","line":"    {"},
{"lineNum":"  386","line":"        template <typename CharEncoding>"},
{"lineNum":"  387","line":"        struct no_case_filter"},
{"lineNum":"  388","line":"        {"},
{"lineNum":"  389","line":"            Char operator()(Char ch) const"},
{"lineNum":"  390","line":"            {"},
{"lineNum":"  391","line":"                return static_cast<Char>(CharEncoding::tolower(ch));"},
{"lineNum":"  392","line":"            }"},
{"lineNum":"  393","line":"        };"},
{"lineNum":"  394","line":""},
{"lineNum":"  395","line":"        typedef has_modifier<Modifiers, tag::char_code_base<tag::no_case> > no_case;"},
{"lineNum":"  396","line":"        typedef reference<symbols<Char, T, Lookup, Filter> > reference_;"},
{"lineNum":"  397","line":"        typedef no_case_filter<"},
{"lineNum":"  398","line":"            typename spirit::detail::get_encoding_with_case<"},
{"lineNum":"  399","line":"                Modifiers"},
{"lineNum":"  400","line":"              , char_encoding::standard"},
{"lineNum":"  401","line":"              , no_case::value>::type>"},
{"lineNum":"  402","line":"        nc_filter;"},
{"lineNum":"  403","line":""},
{"lineNum":"  404","line":"        typedef typename mpl::if_<"},
{"lineNum":"  405","line":"            no_case"},
{"lineNum":"  406","line":"          , symbols<Char, T, Lookup, nc_filter>"},
{"lineNum":"  407","line":"          , reference_>::type"},
{"lineNum":"  408","line":"        result_type;"},
{"lineNum":"  409","line":""},
{"lineNum":"  410","line":"        result_type operator()(reference_ ref, unused_type) const"},
{"lineNum":"  411","line":"        {"},
{"lineNum":"  412","line":"            return result_type(ref.ref.get());"},
{"lineNum":"  413","line":"        }"},
{"lineNum":"  414","line":"    };"},
{"lineNum":"  415","line":"}}}"},
{"lineNum":"  416","line":""},
{"lineNum":"  417","line":"namespace boost { namespace spirit { namespace traits"},
{"lineNum":"  418","line":"{"},
{"lineNum":"  419","line":"    ///////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  420","line":"    template <typename Char, typename T, typename Lookup, typename Filter"},
{"lineNum":"  421","line":"      , typename Attr, typename Context, typename Iterator>"},
{"lineNum":"  422","line":"    struct handles_container<qi::symbols<Char, T, Lookup, Filter>, Attr, Context, Iterator>"},
{"lineNum":"  423","line":"      : traits::is_container<Attr> {};"},
{"lineNum":"  424","line":"}}}"},
{"lineNum":"  425","line":""},
{"lineNum":"  426","line":"#if defined(BOOST_MSVC)"},
{"lineNum":"  427","line":"# pragma warning(pop)"},
{"lineNum":"  428","line":"#endif"},
{"lineNum":"  429","line":""},
{"lineNum":"  430","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 13, "covered" : 0,};
var merged_data = [];
