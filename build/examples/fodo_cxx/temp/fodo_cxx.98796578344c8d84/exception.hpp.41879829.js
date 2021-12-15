var data = {lines:[
{"lineNum":"    1","line":"//Copyright (c) 2006-2009 Emil Dotchevski and Reverge Studios, Inc."},
{"lineNum":"    2","line":""},
{"lineNum":"    3","line":"//Distributed under the Boost Software License, Version 1.0. (See accompanying"},
{"lineNum":"    4","line":"//file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    5","line":""},
{"lineNum":"    6","line":"#ifndef UUID_274DA366004E11DCB1DDFE2E56D89593"},
{"lineNum":"    7","line":"#define UUID_274DA366004E11DCB1DDFE2E56D89593"},
{"lineNum":"    8","line":""},
{"lineNum":"    9","line":"#include <boost/config.hpp>"},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"#ifdef BOOST_EXCEPTION_MINI_BOOST"},
{"lineNum":"   12","line":"#include  <memory>"},
{"lineNum":"   13","line":"namespace boost { namespace exception_detail { using std::shared_ptr; } }"},
{"lineNum":"   14","line":"#else"},
{"lineNum":"   15","line":"namespace boost { template <class T> class shared_ptr; }"},
{"lineNum":"   16","line":"namespace boost { namespace exception_detail { using boost::shared_ptr; } }"},
{"lineNum":"   17","line":"#endif"},
{"lineNum":"   18","line":""},
{"lineNum":"   19","line":"#if defined(__GNUC__) && (__GNUC__*100+__GNUC_MINOR__>301) && !defined(BOOST_EXCEPTION_ENABLE_WARNINGS)"},
{"lineNum":"   20","line":"#pragma GCC system_header"},
{"lineNum":"   21","line":"#endif"},
{"lineNum":"   22","line":"#if defined(_MSC_VER) && !defined(BOOST_EXCEPTION_ENABLE_WARNINGS)"},
{"lineNum":"   23","line":"#pragma warning(push,1)"},
{"lineNum":"   24","line":"#endif"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"namespace"},
{"lineNum":"   27","line":"boost"},
{"lineNum":"   28","line":"    {"},
{"lineNum":"   29","line":"    namespace"},
{"lineNum":"   30","line":"    exception_detail"},
{"lineNum":"   31","line":"        {"},
{"lineNum":"   32","line":"        template <class T>"},
{"lineNum":"   33","line":"        class"},
{"lineNum":"   34","line":"        refcount_ptr"},
{"lineNum":"   35","line":"            {"},
{"lineNum":"   36","line":"            public:"},
{"lineNum":"   37","line":""},
{"lineNum":"   38","line":"            refcount_ptr():"},
{"lineNum":"   39","line":"                px_(0)"},
{"lineNum":"   40","line":"                {"},
{"lineNum":"   41","line":"                }"},
{"lineNum":"   42","line":""},
{"lineNum":"   43","line":"            ~refcount_ptr()"},
{"lineNum":"   44","line":"                {"},
{"lineNum":"   45","line":"                release();","class":"lineNoCov","hits":"0","possible_hits":"109",},
{"lineNum":"   46","line":"                }"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"            refcount_ptr( refcount_ptr const & x ):"},
{"lineNum":"   49","line":"                px_(x.px_)","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   50","line":"                {"},
{"lineNum":"   51","line":"                add_ref();"},
{"lineNum":"   52","line":"                }"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":"            refcount_ptr &"},
{"lineNum":"   55","line":"            operator=( refcount_ptr const & x )"},
{"lineNum":"   56","line":"                {"},
{"lineNum":"   57","line":"                adopt(x.px_);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   58","line":"                return *this;"},
{"lineNum":"   59","line":"                }"},
{"lineNum":"   60","line":""},
{"lineNum":"   61","line":"            void"},
{"lineNum":"   62","line":"            adopt( T * px )"},
{"lineNum":"   63","line":"                {"},
{"lineNum":"   64","line":"                release();"},
{"lineNum":"   65","line":"                px_=px;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   66","line":"                add_ref();"},
{"lineNum":"   67","line":"                }"},
{"lineNum":"   68","line":""},
{"lineNum":"   69","line":"            T *"},
{"lineNum":"   70","line":"            get() const"},
{"lineNum":"   71","line":"                {"},
{"lineNum":"   72","line":"                return px_;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   73","line":"                }"},
{"lineNum":"   74","line":""},
{"lineNum":"   75","line":"            private:"},
{"lineNum":"   76","line":""},
{"lineNum":"   77","line":"            T * px_;"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":"            void"},
{"lineNum":"   80","line":"            add_ref()"},
{"lineNum":"   81","line":"                {"},
{"lineNum":"   82","line":"                if( px_ )","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   83","line":"                    px_->add_ref();","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"   84","line":"                }"},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"            void"},
{"lineNum":"   87","line":"            release()"},
{"lineNum":"   88","line":"                {"},
{"lineNum":"   89","line":"                if( px_ && px_->release() )","class":"lineNoCov","hits":"0","possible_hits":"59",},
{"lineNum":"   90","line":"                    px_=0;","class":"lineNoCov","hits":"0","possible_hits":"54",},
{"lineNum":"   91","line":"                }"},
{"lineNum":"   92","line":"            };"},
{"lineNum":"   93","line":"        }"},
{"lineNum":"   94","line":""},
{"lineNum":"   95","line":"    ////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   96","line":""},
{"lineNum":"   97","line":"    template <class Tag,class T>"},
{"lineNum":"   98","line":"    class error_info;"},
{"lineNum":"   99","line":""},
{"lineNum":"  100","line":"    typedef error_info<struct throw_function_,char const *> throw_function;"},
{"lineNum":"  101","line":"    typedef error_info<struct throw_file_,char const *> throw_file;"},
{"lineNum":"  102","line":"    typedef error_info<struct throw_line_,int> throw_line;"},
{"lineNum":"  103","line":""},
{"lineNum":"  104","line":"    template <>"},
{"lineNum":"  105","line":"    class"},
{"lineNum":"  106","line":"    error_info<throw_function_,char const *>"},
{"lineNum":"  107","line":"        {"},
{"lineNum":"  108","line":"        public:"},
{"lineNum":"  109","line":"        typedef char const * value_type;"},
{"lineNum":"  110","line":"        value_type v_;"},
{"lineNum":"  111","line":"        explicit"},
{"lineNum":"  112","line":"        error_info( value_type v ):"},
{"lineNum":"  113","line":"            v_(v)"},
{"lineNum":"  114","line":"            {"},
{"lineNum":"  115","line":"            }"},
{"lineNum":"  116","line":"        };"},
{"lineNum":"  117","line":""},
{"lineNum":"  118","line":"    template <>"},
{"lineNum":"  119","line":"    class"},
{"lineNum":"  120","line":"    error_info<throw_file_,char const *>"},
{"lineNum":"  121","line":"        {"},
{"lineNum":"  122","line":"        public:"},
{"lineNum":"  123","line":"        typedef char const * value_type;"},
{"lineNum":"  124","line":"        value_type v_;"},
{"lineNum":"  125","line":"        explicit"},
{"lineNum":"  126","line":"        error_info( value_type v ):"},
{"lineNum":"  127","line":"            v_(v)"},
{"lineNum":"  128","line":"            {"},
{"lineNum":"  129","line":"            }"},
{"lineNum":"  130","line":"        };"},
{"lineNum":"  131","line":""},
{"lineNum":"  132","line":"    template <>"},
{"lineNum":"  133","line":"    class"},
{"lineNum":"  134","line":"    error_info<throw_line_,int>"},
{"lineNum":"  135","line":"        {"},
{"lineNum":"  136","line":"        public:"},
{"lineNum":"  137","line":"        typedef int value_type;"},
{"lineNum":"  138","line":"        value_type v_;"},
{"lineNum":"  139","line":"        explicit"},
{"lineNum":"  140","line":"        error_info( value_type v ):"},
{"lineNum":"  141","line":"            v_(v)"},
{"lineNum":"  142","line":"            {"},
{"lineNum":"  143","line":"            }"},
{"lineNum":"  144","line":"        };"},
{"lineNum":"  145","line":""},
{"lineNum":"  146","line":"    class"},
{"lineNum":"  147","line":"    BOOST_SYMBOL_VISIBLE"},
{"lineNum":"  148","line":"    exception;"},
{"lineNum":"  149","line":""},
{"lineNum":"  150","line":"    namespace"},
{"lineNum":"  151","line":"    exception_detail"},
{"lineNum":"  152","line":"        {"},
{"lineNum":"  153","line":"        class error_info_base;"},
{"lineNum":"  154","line":"        struct type_info_;"},
{"lineNum":"  155","line":""},
{"lineNum":"  156","line":"        struct"},
{"lineNum":"  157","line":"        error_info_container"},
{"lineNum":"  158","line":"            {"},
{"lineNum":"  159","line":"            virtual char const * diagnostic_information( char const * ) const = 0;"},
{"lineNum":"  160","line":"            virtual shared_ptr<error_info_base> get( type_info_ const & ) const = 0;"},
{"lineNum":"  161","line":"            virtual void set( shared_ptr<error_info_base> const &, type_info_ const & ) = 0;"},
{"lineNum":"  162","line":"            virtual void add_ref() const = 0;"},
{"lineNum":"  163","line":"            virtual bool release() const = 0;"},
{"lineNum":"  164","line":"            virtual refcount_ptr<exception_detail::error_info_container> clone() const = 0;"},
{"lineNum":"  165","line":""},
{"lineNum":"  166","line":"            protected:"},
{"lineNum":"  167","line":""},
{"lineNum":"  168","line":"            ~error_info_container() throw()"},
{"lineNum":"  169","line":"                {"},
{"lineNum":"  170","line":"                }"},
{"lineNum":"  171","line":"            };"},
{"lineNum":"  172","line":""},
{"lineNum":"  173","line":"        template <class>"},
{"lineNum":"  174","line":"        struct get_info;"},
{"lineNum":"  175","line":""},
{"lineNum":"  176","line":"        template <>"},
{"lineNum":"  177","line":"        struct get_info<throw_function>;"},
{"lineNum":"  178","line":""},
{"lineNum":"  179","line":"        template <>"},
{"lineNum":"  180","line":"        struct get_info<throw_file>;"},
{"lineNum":"  181","line":""},
{"lineNum":"  182","line":"        template <>"},
{"lineNum":"  183","line":"        struct get_info<throw_line>;"},
{"lineNum":"  184","line":""},
{"lineNum":"  185","line":"        template <class>"},
{"lineNum":"  186","line":"        struct set_info_rv;"},
{"lineNum":"  187","line":""},
{"lineNum":"  188","line":"        template <>"},
{"lineNum":"  189","line":"        struct set_info_rv<throw_function>;"},
{"lineNum":"  190","line":""},
{"lineNum":"  191","line":"        template <>"},
{"lineNum":"  192","line":"        struct set_info_rv<throw_file>;"},
{"lineNum":"  193","line":""},
{"lineNum":"  194","line":"        template <>"},
{"lineNum":"  195","line":"        struct set_info_rv<throw_line>;"},
{"lineNum":"  196","line":""},
{"lineNum":"  197","line":"        char const * get_diagnostic_information( exception const &, char const * );"},
{"lineNum":"  198","line":""},
{"lineNum":"  199","line":"        void copy_boost_exception( exception *, exception const * );"},
{"lineNum":"  200","line":""},
{"lineNum":"  201","line":"        template <class E,class Tag,class T>"},
{"lineNum":"  202","line":"        E const & set_info( E const &, error_info<Tag,T> const & );"},
{"lineNum":"  203","line":""},
{"lineNum":"  204","line":"        template <class E>"},
{"lineNum":"  205","line":"        E const & set_info( E const &, throw_function const & );"},
{"lineNum":"  206","line":""},
{"lineNum":"  207","line":"        template <class E>"},
{"lineNum":"  208","line":"        E const & set_info( E const &, throw_file const & );"},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"        template <class E>"},
{"lineNum":"  211","line":"        E const & set_info( E const &, throw_line const & );"},
{"lineNum":"  212","line":"        }"},
{"lineNum":"  213","line":""},
{"lineNum":"  214","line":"    class"},
{"lineNum":"  215","line":"    BOOST_SYMBOL_VISIBLE"},
{"lineNum":"  216","line":"    exception","class":"lineNoCov","hits":"0","possible_hits":"22",},
{"lineNum":"  217","line":"        {"},
{"lineNum":"  218","line":"        //<N3757>"},
{"lineNum":"  219","line":"        public:"},
{"lineNum":"  220","line":"        template <class Tag> void set( typename Tag::type const & );"},
{"lineNum":"  221","line":"        template <class Tag> typename Tag::type const * get() const;"},
{"lineNum":"  222","line":"        //</N3757>"},
{"lineNum":"  223","line":""},
{"lineNum":"  224","line":"        protected:"},
{"lineNum":"  225","line":""},
{"lineNum":"  226","line":"        exception():"},
{"lineNum":"  227","line":"            throw_function_(0),","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  228","line":"            throw_file_(0),"},
{"lineNum":"  229","line":"            throw_line_(-1)","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  230","line":"            {","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  231","line":"            }"},
{"lineNum":"  232","line":""},
{"lineNum":"  233","line":"#ifdef __HP_aCC"},
{"lineNum":"  234","line":"        //On HP aCC, this protected copy constructor prevents throwing boost::exception."},
{"lineNum":"  235","line":"        //On all other platforms, the same effect is achieved by the pure virtual destructor."},
{"lineNum":"  236","line":"        exception( exception const & x ) throw():"},
{"lineNum":"  237","line":"            data_(x.data_),"},
{"lineNum":"  238","line":"            throw_function_(x.throw_function_),"},
{"lineNum":"  239","line":"            throw_file_(x.throw_file_),"},
{"lineNum":"  240","line":"            throw_line_(x.throw_line_)"},
{"lineNum":"  241","line":"            {"},
{"lineNum":"  242","line":"            }"},
{"lineNum":"  243","line":"#endif"},
{"lineNum":"  244","line":""},
{"lineNum":"  245","line":"        virtual ~exception() throw()"},
{"lineNum":"  246","line":"#ifndef __HP_aCC"},
{"lineNum":"  247","line":"            = 0 //Workaround for HP aCC, =0 incorrectly leads to link errors."},
{"lineNum":"  248","line":"#endif"},
{"lineNum":"  249","line":"            ;"},
{"lineNum":"  250","line":""},
{"lineNum":"  251","line":"#if (defined(__MWERKS__) && __MWERKS__<=0x3207) || (defined(_MSC_VER) && _MSC_VER<=1310)"},
{"lineNum":"  252","line":"        public:"},
{"lineNum":"  253","line":"#else"},
{"lineNum":"  254","line":"        private:"},
{"lineNum":"  255","line":""},
{"lineNum":"  256","line":"        template <class E>"},
{"lineNum":"  257","line":"        friend E const & exception_detail::set_info( E const &, throw_function const & );"},
{"lineNum":"  258","line":""},
{"lineNum":"  259","line":"        template <class E>"},
{"lineNum":"  260","line":"        friend E const & exception_detail::set_info( E const &, throw_file const & );"},
{"lineNum":"  261","line":""},
{"lineNum":"  262","line":"        template <class E>"},
{"lineNum":"  263","line":"        friend E const & exception_detail::set_info( E const &, throw_line const & );"},
{"lineNum":"  264","line":""},
{"lineNum":"  265","line":"        template <class E,class Tag,class T>"},
{"lineNum":"  266","line":"        friend E const & exception_detail::set_info( E const &, error_info<Tag,T> const & );"},
{"lineNum":"  267","line":""},
{"lineNum":"  268","line":"        friend char const * exception_detail::get_diagnostic_information( exception const &, char const * );"},
{"lineNum":"  269","line":""},
{"lineNum":"  270","line":"        template <class>"},
{"lineNum":"  271","line":"        friend struct exception_detail::get_info;"},
{"lineNum":"  272","line":"        friend struct exception_detail::get_info<throw_function>;"},
{"lineNum":"  273","line":"        friend struct exception_detail::get_info<throw_file>;"},
{"lineNum":"  274","line":"        friend struct exception_detail::get_info<throw_line>;"},
{"lineNum":"  275","line":"        template <class>"},
{"lineNum":"  276","line":"        friend struct exception_detail::set_info_rv;"},
{"lineNum":"  277","line":"        friend struct exception_detail::set_info_rv<throw_function>;"},
{"lineNum":"  278","line":"        friend struct exception_detail::set_info_rv<throw_file>;"},
{"lineNum":"  279","line":"        friend struct exception_detail::set_info_rv<throw_line>;"},
{"lineNum":"  280","line":"        friend void exception_detail::copy_boost_exception( exception *, exception const * );"},
{"lineNum":"  281","line":"#endif"},
{"lineNum":"  282","line":"        mutable exception_detail::refcount_ptr<exception_detail::error_info_container> data_;"},
{"lineNum":"  283","line":"        mutable char const * throw_function_;"},
{"lineNum":"  284","line":"        mutable char const * throw_file_;"},
{"lineNum":"  285","line":"        mutable int throw_line_;"},
{"lineNum":"  286","line":"        };"},
{"lineNum":"  287","line":""},
{"lineNum":"  288","line":"    inline"},
{"lineNum":"  289","line":"    exception::"},
{"lineNum":"  290","line":"    ~exception() throw()"},
{"lineNum":"  291","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"54",},
{"lineNum":"  292","line":"        }"},
{"lineNum":"  293","line":""},
{"lineNum":"  294","line":"    namespace"},
{"lineNum":"  295","line":"    exception_detail"},
{"lineNum":"  296","line":"        {"},
{"lineNum":"  297","line":"        template <class E>"},
{"lineNum":"  298","line":"        E const &"},
{"lineNum":"  299","line":"        set_info( E const & x, throw_function const & y )"},
{"lineNum":"  300","line":"            {"},
{"lineNum":"  301","line":"            x.throw_function_=y.v_;"},
{"lineNum":"  302","line":"            return x;"},
{"lineNum":"  303","line":"            }"},
{"lineNum":"  304","line":""},
{"lineNum":"  305","line":"        template <class E>"},
{"lineNum":"  306","line":"        E const &"},
{"lineNum":"  307","line":"        set_info( E const & x, throw_file const & y )"},
{"lineNum":"  308","line":"            {"},
{"lineNum":"  309","line":"            x.throw_file_=y.v_;"},
{"lineNum":"  310","line":"            return x;"},
{"lineNum":"  311","line":"            }"},
{"lineNum":"  312","line":""},
{"lineNum":"  313","line":"        template <class E>"},
{"lineNum":"  314","line":"        E const &"},
{"lineNum":"  315","line":"        set_info( E const & x, throw_line const & y )"},
{"lineNum":"  316","line":"            {"},
{"lineNum":"  317","line":"            x.throw_line_=y.v_;"},
{"lineNum":"  318","line":"            return x;"},
{"lineNum":"  319","line":"            }"},
{"lineNum":"  320","line":"        }"},
{"lineNum":"  321","line":""},
{"lineNum":"  322","line":"    ////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  323","line":""},
{"lineNum":"  324","line":"    namespace"},
{"lineNum":"  325","line":"    exception_detail"},
{"lineNum":"  326","line":"        {"},
{"lineNum":"  327","line":"        template <class T>"},
{"lineNum":"  328","line":"        struct"},
{"lineNum":"  329","line":"        BOOST_SYMBOL_VISIBLE"},
{"lineNum":"  330","line":"        error_info_injector:","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  331","line":"            public T,"},
{"lineNum":"  332","line":"            public exception"},
{"lineNum":"  333","line":"            {"},
{"lineNum":"  334","line":"            explicit"},
{"lineNum":"  335","line":"            error_info_injector( T const & x ):"},
{"lineNum":"  336","line":"                T(x)"},
{"lineNum":"  337","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  338","line":"                }"},
{"lineNum":"  339","line":""},
{"lineNum":"  340","line":"            ~error_info_injector() throw()"},
{"lineNum":"  341","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"  342","line":"                }","class":"lineNoCov","hits":"0","possible_hits":"51",},
{"lineNum":"  343","line":"            };"},
{"lineNum":"  344","line":""},
{"lineNum":"  345","line":"        struct large_size { char c[256]; };"},
{"lineNum":"  346","line":"        large_size dispatch_boost_exception( exception const * );"},
{"lineNum":"  347","line":""},
{"lineNum":"  348","line":"        struct small_size { };"},
{"lineNum":"  349","line":"        small_size dispatch_boost_exception( void const * );"},
{"lineNum":"  350","line":""},
{"lineNum":"  351","line":"        template <class,int>"},
{"lineNum":"  352","line":"        struct enable_error_info_helper;"},
{"lineNum":"  353","line":""},
{"lineNum":"  354","line":"        template <class T>"},
{"lineNum":"  355","line":"        struct"},
{"lineNum":"  356","line":"        enable_error_info_helper<T,sizeof(large_size)>"},
{"lineNum":"  357","line":"            {"},
{"lineNum":"  358","line":"            typedef T type;"},
{"lineNum":"  359","line":"            };"},
{"lineNum":"  360","line":""},
{"lineNum":"  361","line":"        template <class T>"},
{"lineNum":"  362","line":"        struct"},
{"lineNum":"  363","line":"        enable_error_info_helper<T,sizeof(small_size)>"},
{"lineNum":"  364","line":"            {"},
{"lineNum":"  365","line":"            typedef error_info_injector<T> type;"},
{"lineNum":"  366","line":"            };"},
{"lineNum":"  367","line":""},
{"lineNum":"  368","line":"        template <class T>"},
{"lineNum":"  369","line":"        struct"},
{"lineNum":"  370","line":"        enable_error_info_return_type"},
{"lineNum":"  371","line":"            {"},
{"lineNum":"  372","line":"            typedef typename enable_error_info_helper<T,sizeof(exception_detail::dispatch_boost_exception(static_cast<T *>(0)))>::type type;"},
{"lineNum":"  373","line":"            };"},
{"lineNum":"  374","line":"        }"},
{"lineNum":"  375","line":""},
{"lineNum":"  376","line":"    template <class T>"},
{"lineNum":"  377","line":"    inline"},
{"lineNum":"  378","line":"    typename"},
{"lineNum":"  379","line":"    exception_detail::enable_error_info_return_type<T>::type"},
{"lineNum":"  380","line":"    enable_error_info( T const & x )"},
{"lineNum":"  381","line":"        {"},
{"lineNum":"  382","line":"        typedef typename exception_detail::enable_error_info_return_type<T>::type rt;"},
{"lineNum":"  383","line":"        return rt(x);"},
{"lineNum":"  384","line":"        }"},
{"lineNum":"  385","line":""},
{"lineNum":"  386","line":"    ////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  387","line":""},
{"lineNum":"  388","line":"    namespace"},
{"lineNum":"  389","line":"    exception_detail"},
{"lineNum":"  390","line":"        {"},
{"lineNum":"  391","line":"        class"},
{"lineNum":"  392","line":"        BOOST_SYMBOL_VISIBLE"},
{"lineNum":"  393","line":"        clone_base","class":"lineNoCov","hits":"0","possible_hits":"10",},
{"lineNum":"  394","line":"            {"},
{"lineNum":"  395","line":"            public:"},
{"lineNum":"  396","line":""},
{"lineNum":"  397","line":"            virtual clone_base const * clone() const = 0;"},
{"lineNum":"  398","line":"            virtual void rethrow() const = 0;"},
{"lineNum":"  399","line":""},
{"lineNum":"  400","line":"            virtual"},
{"lineNum":"  401","line":"            ~clone_base() throw()"},
{"lineNum":"  402","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  403","line":"                }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  404","line":"            };"},
{"lineNum":"  405","line":""},
{"lineNum":"  406","line":"        inline"},
{"lineNum":"  407","line":"        void"},
{"lineNum":"  408","line":"        copy_boost_exception( exception * a, exception const * b )"},
{"lineNum":"  409","line":"            {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  410","line":"            refcount_ptr<error_info_container> data;"},
{"lineNum":"  411","line":"            if( error_info_container * d=b->data_.get() )","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  412","line":"                data = d->clone();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  413","line":"            a->throw_file_ = b->throw_file_;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  414","line":"            a->throw_line_ = b->throw_line_;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  415","line":"            a->throw_function_ = b->throw_function_;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  416","line":"            a->data_ = data;"},
{"lineNum":"  417","line":"            }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  418","line":""},
{"lineNum":"  419","line":"        inline"},
{"lineNum":"  420","line":"        void"},
{"lineNum":"  421","line":"        copy_boost_exception( void *, void const * )"},
{"lineNum":"  422","line":"            {"},
{"lineNum":"  423","line":"            }"},
{"lineNum":"  424","line":""},
{"lineNum":"  425","line":"        template <class T>"},
{"lineNum":"  426","line":"        class"},
{"lineNum":"  427","line":"        BOOST_SYMBOL_VISIBLE"},
{"lineNum":"  428","line":"        clone_impl:","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  429","line":"            public T,"},
{"lineNum":"  430","line":"            public virtual clone_base"},
{"lineNum":"  431","line":"            {"},
{"lineNum":"  432","line":"            struct clone_tag { };"},
{"lineNum":"  433","line":"            clone_impl( clone_impl const & x, clone_tag ):"},
{"lineNum":"  434","line":"                T(x)"},
{"lineNum":"  435","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  436","line":"                copy_boost_exception(this,&x);","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  437","line":"                }","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  438","line":""},
{"lineNum":"  439","line":"            public:"},
{"lineNum":"  440","line":""},
{"lineNum":"  441","line":"            explicit"},
{"lineNum":"  442","line":"            clone_impl( T const & x ):"},
{"lineNum":"  443","line":"                T(x)"},
{"lineNum":"  444","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  445","line":"                copy_boost_exception(this,&x);","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  446","line":"                }","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  447","line":""},
{"lineNum":"  448","line":"            ~clone_impl() throw()"},
{"lineNum":"  449","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"15",},
{"lineNum":"  450","line":"                }"},
{"lineNum":"  451","line":""},
{"lineNum":"  452","line":"            private:"},
{"lineNum":"  453","line":""},
{"lineNum":"  454","line":"            clone_base const *"},
{"lineNum":"  455","line":"            clone() const"},
{"lineNum":"  456","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  457","line":"                return new clone_impl(*this,clone_tag());","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"  458","line":"                }"},
{"lineNum":"  459","line":""},
{"lineNum":"  460","line":"            void"},
{"lineNum":"  461","line":"            rethrow() const"},
{"lineNum":"  462","line":"                {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  463","line":"                throw*this;","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  464","line":"                }"},
{"lineNum":"  465","line":"            };"},
{"lineNum":"  466","line":"        }"},
{"lineNum":"  467","line":""},
{"lineNum":"  468","line":"    template <class T>"},
{"lineNum":"  469","line":"    inline"},
{"lineNum":"  470","line":"    exception_detail::clone_impl<T>"},
{"lineNum":"  471","line":"    enable_current_exception( T const & x )"},
{"lineNum":"  472","line":"        {"},
{"lineNum":"  473","line":"        return exception_detail::clone_impl<T>(x);"},
{"lineNum":"  474","line":"        }"},
{"lineNum":"  475","line":""},
{"lineNum":"  476","line":"    template <class T>"},
{"lineNum":"  477","line":"    struct"},
{"lineNum":"  478","line":"    BOOST_SYMBOL_VISIBLE"},
{"lineNum":"  479","line":"    wrapexcept:"},
{"lineNum":"  480","line":"        public exception_detail::clone_impl<typename exception_detail::enable_error_info_return_type<T>::type>"},
{"lineNum":"  481","line":"        {"},
{"lineNum":"  482","line":"        typedef exception_detail::clone_impl<typename exception_detail::enable_error_info_return_type<T>::type> base_type;"},
{"lineNum":"  483","line":"        public:"},
{"lineNum":"  484","line":"        explicit"},
{"lineNum":"  485","line":"        wrapexcept( typename exception_detail::enable_error_info_return_type<T>::type const & x ):"},
{"lineNum":"  486","line":"            base_type( x )"},
{"lineNum":"  487","line":"            {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  488","line":"            }"},
{"lineNum":"  489","line":""},
{"lineNum":"  490","line":"        ~wrapexcept() throw()"},
{"lineNum":"  491","line":"            {","class":"lineNoCov","hits":"0","possible_hits":"15",},
{"lineNum":"  492","line":"            }"},
{"lineNum":"  493","line":"        };"},
{"lineNum":"  494","line":""},
{"lineNum":"  495","line":"    namespace"},
{"lineNum":"  496","line":"    exception_detail"},
{"lineNum":"  497","line":"        {"},
{"lineNum":"  498","line":"        template <class T>"},
{"lineNum":"  499","line":"        struct"},
{"lineNum":"  500","line":"        remove_error_info_injector"},
{"lineNum":"  501","line":"            {"},
{"lineNum":"  502","line":"            typedef T type;"},
{"lineNum":"  503","line":"            };"},
{"lineNum":"  504","line":""},
{"lineNum":"  505","line":"        template <class T>"},
{"lineNum":"  506","line":"        struct"},
{"lineNum":"  507","line":"        remove_error_info_injector< error_info_injector<T> >"},
{"lineNum":"  508","line":"            {"},
{"lineNum":"  509","line":"            typedef T type;"},
{"lineNum":"  510","line":"            };"},
{"lineNum":"  511","line":""},
{"lineNum":"  512","line":"        template <class T>"},
{"lineNum":"  513","line":"        inline"},
{"lineNum":"  514","line":"        wrapexcept<typename remove_error_info_injector<T>::type>"},
{"lineNum":"  515","line":"        enable_both( T const & x )"},
{"lineNum":"  516","line":"            {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  517","line":"            return wrapexcept<typename remove_error_info_injector<T>::type>( enable_error_info( x ) );","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  518","line":"            }"},
{"lineNum":"  519","line":"        }"},
{"lineNum":"  520","line":"    }"},
{"lineNum":"  521","line":""},
{"lineNum":"  522","line":"#if defined(_MSC_VER) && !defined(BOOST_EXCEPTION_ENABLE_WARNINGS)"},
{"lineNum":"  523","line":"#pragma warning(pop)"},
{"lineNum":"  524","line":"#endif"},
{"lineNum":"  525","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 44, "covered" : 0,};
var merged_data = [];
