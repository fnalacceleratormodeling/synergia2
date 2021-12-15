var data = {lines:[
{"lineNum":"    1","line":"/*! \\file access.hpp"},
{"lineNum":"    2","line":"    \\brief Access control and default construction */"},
{"lineNum":"    3","line":"/*"},
{"lineNum":"    4","line":"  Copyright (c) 2014, Randolph Voorhies, Shane Grant"},
{"lineNum":"    5","line":"  All rights reserved."},
{"lineNum":"    6","line":""},
{"lineNum":"    7","line":"  Redistribution and use in source and binary forms, with or without"},
{"lineNum":"    8","line":"  modification, are permitted provided that the following conditions are met:"},
{"lineNum":"    9","line":"      * Redistributions of source code must retain the above copyright"},
{"lineNum":"   10","line":"        notice, this list of conditions and the following disclaimer."},
{"lineNum":"   11","line":"      * Redistributions in binary form must reproduce the above copyright"},
{"lineNum":"   12","line":"        notice, this list of conditions and the following disclaimer in the"},
{"lineNum":"   13","line":"        documentation and/or other materials provided with the distribution."},
{"lineNum":"   14","line":"      * Neither the name of cereal nor the"},
{"lineNum":"   15","line":"        names of its contributors may be used to endorse or promote products"},
{"lineNum":"   16","line":"        derived from this software without specific prior written permission."},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND"},
{"lineNum":"   19","line":"  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED"},
{"lineNum":"   20","line":"  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE"},
{"lineNum":"   21","line":"  DISCLAIMED. IN NO EVENT SHALL RANDOLPH VOORHIES OR SHANE GRANT BE LIABLE FOR ANY"},
{"lineNum":"   22","line":"  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES"},
{"lineNum":"   23","line":"  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;"},
{"lineNum":"   24","line":"  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND"},
{"lineNum":"   25","line":"  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT"},
{"lineNum":"   26","line":"  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"},
{"lineNum":"   27","line":"  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."},
{"lineNum":"   28","line":"*/"},
{"lineNum":"   29","line":"#ifndef CEREAL_ACCESS_HPP_"},
{"lineNum":"   30","line":"#define CEREAL_ACCESS_HPP_"},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":"#include <type_traits>"},
{"lineNum":"   33","line":"#include <iostream>"},
{"lineNum":"   34","line":"#include <cstdint>"},
{"lineNum":"   35","line":"#include <functional>"},
{"lineNum":"   36","line":""},
{"lineNum":"   37","line":"#include \"cereal/macros.hpp\""},
{"lineNum":"   38","line":"#include \"cereal/specialize.hpp\""},
{"lineNum":"   39","line":"#include \"cereal/details/helpers.hpp\""},
{"lineNum":"   40","line":""},
{"lineNum":"   41","line":"namespace cereal"},
{"lineNum":"   42","line":"{"},
{"lineNum":"   43","line":"  // ######################################################################"},
{"lineNum":"   44","line":"  //! A class that allows cereal to load smart pointers to types that have no default constructor"},
{"lineNum":"   45","line":"  /*! If your class does not have a default constructor, cereal will not be able"},
{"lineNum":"   46","line":"      to load any smart pointers to it unless you overload LoadAndConstruct"},
{"lineNum":"   47","line":"      for your class, and provide an appropriate load_and_construct method.  You can also"},
{"lineNum":"   48","line":"      choose to define a member static function instead of specializing this class."},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"      The specialization of LoadAndConstruct must be placed within the cereal namespace:"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"      @code{.cpp}"},
{"lineNum":"   53","line":"      struct MyType"},
{"lineNum":"   54","line":"      {"},
{"lineNum":"   55","line":"        MyType( int x ); // note: no default ctor"},
{"lineNum":"   56","line":"        int myX;"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"        // Define a serialize or load/save pair as you normally would"},
{"lineNum":"   59","line":"        template <class Archive>"},
{"lineNum":"   60","line":"        void serialize( Archive & ar )"},
{"lineNum":"   61","line":"        {"},
{"lineNum":"   62","line":"          ar( myX );"},
{"lineNum":"   63","line":"        }"},
{"lineNum":"   64","line":"      };"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"      // Provide a specialization for LoadAndConstruct for your type"},
{"lineNum":"   67","line":"      namespace cereal"},
{"lineNum":"   68","line":"      {"},
{"lineNum":"   69","line":"        template <> struct LoadAndConstruct<MyType>"},
{"lineNum":"   70","line":"        {"},
{"lineNum":"   71","line":"          // load_and_construct will be passed the archive that you will be loading"},
{"lineNum":"   72","line":"          // from as well as a construct object which you can use as if it were the"},
{"lineNum":"   73","line":"          // constructor for your type.  cereal will handle all memory management for you."},
{"lineNum":"   74","line":"          template <class Archive>"},
{"lineNum":"   75","line":"          static void load_and_construct( Archive & ar, cereal::construct<MyType> & construct )"},
{"lineNum":"   76","line":"          {"},
{"lineNum":"   77","line":"            int x;"},
{"lineNum":"   78","line":"            ar( x );"},
{"lineNum":"   79","line":"            construct( x );"},
{"lineNum":"   80","line":"          }"},
{"lineNum":"   81","line":""},
{"lineNum":"   82","line":"          // if you require versioning, simply add a const std::uint32_t as the final parameter, e.g.:"},
{"lineNum":"   83","line":"          // load_and_construct( Archive & ar, cereal::construct<MyType> & construct, std::uint32_t const version )"},
{"lineNum":"   84","line":"        };"},
{"lineNum":"   85","line":"      } // end namespace cereal"},
{"lineNum":"   86","line":"      @endcode"},
{"lineNum":"   87","line":""},
{"lineNum":"   88","line":"      Please note that just as in using external serialization functions, you cannot get"},
{"lineNum":"   89","line":"      access to non-public members of your class by befriending cereal::access.  If you"},
{"lineNum":"   90","line":"      have the ability to modify the class you wish to serialize, it is recommended that you"},
{"lineNum":"   91","line":"      use member serialize functions and a static member load_and_construct function."},
{"lineNum":"   92","line":""},
{"lineNum":"   93","line":"      load_and_construct functions, regardless of whether they are static members of your class or"},
{"lineNum":"   94","line":"      whether you create one in the LoadAndConstruct specialization, have the following signature:"},
{"lineNum":"   95","line":""},
{"lineNum":"   96","line":"      @code{.cpp}"},
{"lineNum":"   97","line":"      // generally Archive will be templated, but it can be specific if desired"},
{"lineNum":"   98","line":"      template <class Archive>"},
{"lineNum":"   99","line":"      static void load_and_construct( Archive & ar, cereal::construct<MyType> & construct );"},
{"lineNum":"  100","line":"      // with an optional last parameter specifying the version: const std::uint32_t version"},
{"lineNum":"  101","line":"      @endcode"},
{"lineNum":"  102","line":""},
{"lineNum":"  103","line":"      Versioning behaves the same way as it does for standard serialization functions."},
{"lineNum":"  104","line":""},
{"lineNum":"  105","line":"      @tparam T The type to specialize for"},
{"lineNum":"  106","line":"      @ingroup Access */"},
{"lineNum":"  107","line":"  template <class T>"},
{"lineNum":"  108","line":"  struct LoadAndConstruct"},
{"lineNum":"  109","line":"  { };"},
{"lineNum":"  110","line":""},
{"lineNum":"  111","line":"  // forward decl for construct"},
{"lineNum":"  112","line":"  //! @cond PRIVATE_NEVERDEFINED"},
{"lineNum":"  113","line":"  namespace memory_detail{ template <class Ar, class T> struct LoadAndConstructLoadWrapper; }"},
{"lineNum":"  114","line":"  namespace boost_variant_detail{ template <class Ar, class T> struct LoadAndConstructLoadWrapper; }"},
{"lineNum":"  115","line":"  //! @endcond"},
{"lineNum":"  116","line":""},
{"lineNum":"  117","line":"  //! Used to construct types with no default constructor"},
{"lineNum":"  118","line":"  /*! When serializing a type that has no default constructor, cereal"},
{"lineNum":"  119","line":"      will attempt to call either the class static function load_and_construct"},
{"lineNum":"  120","line":"      or the appropriate template specialization of LoadAndConstruct.  cereal"},
{"lineNum":"  121","line":"      will pass that function a reference to the archive as well as a reference"},
{"lineNum":"  122","line":"      to a construct object which should be used to perform the allocation once"},
{"lineNum":"  123","line":"      data has been appropriately loaded."},
{"lineNum":"  124","line":""},
{"lineNum":"  125","line":"      @code{.cpp}"},
{"lineNum":"  126","line":"      struct MyType"},
{"lineNum":"  127","line":"      {"},
{"lineNum":"  128","line":"        // note the lack of default constructor"},
{"lineNum":"  129","line":"        MyType( int xx, int yy );"},
{"lineNum":"  130","line":""},
{"lineNum":"  131","line":"        int x, y;"},
{"lineNum":"  132","line":"        double notInConstructor;"},
{"lineNum":"  133","line":""},
{"lineNum":"  134","line":"        template <class Archive>"},
{"lineNum":"  135","line":"        void serialize( Archive & ar )"},
{"lineNum":"  136","line":"        {"},
{"lineNum":"  137","line":"          ar( x, y );"},
{"lineNum":"  138","line":"          ar( notInConstructor );"},
{"lineNum":"  139","line":"        }"},
{"lineNum":"  140","line":""},
{"lineNum":"  141","line":"        template <class Archive>"},
{"lineNum":"  142","line":"        static void load_and_construct( Archive & ar, cereal::construct<MyType> & construct )"},
{"lineNum":"  143","line":"        {"},
{"lineNum":"  144","line":"          int x, y;"},
{"lineNum":"  145","line":"          ar( x, y );"},
{"lineNum":"  146","line":""},
{"lineNum":"  147","line":"          // use construct object to initialize with loaded data"},
{"lineNum":"  148","line":"          construct( x, y );"},
{"lineNum":"  149","line":""},
{"lineNum":"  150","line":"          // access to member variables and functions via -> operator"},
{"lineNum":"  151","line":"          ar( construct->notInConstructor );"},
{"lineNum":"  152","line":""},
{"lineNum":"  153","line":"          // could also do the above section by:"},
{"lineNum":"  154","line":"          double z;"},
{"lineNum":"  155","line":"          ar( z );"},
{"lineNum":"  156","line":"          construct->notInConstructor = z;"},
{"lineNum":"  157","line":"        }"},
{"lineNum":"  158","line":"      };"},
{"lineNum":"  159","line":"      @endcode"},
{"lineNum":"  160","line":""},
{"lineNum":"  161","line":"      @tparam T The class type being serialized"},
{"lineNum":"  162","line":"      */"},
{"lineNum":"  163","line":"  template <class T>"},
{"lineNum":"  164","line":"  class construct"},
{"lineNum":"  165","line":"  {"},
{"lineNum":"  166","line":"    public:"},
{"lineNum":"  167","line":"      //! Construct and initialize the type T with the given arguments"},
{"lineNum":"  168","line":"      /*! This will forward all arguments to the underlying type T,"},
{"lineNum":"  169","line":"          calling an appropriate constructor."},
{"lineNum":"  170","line":""},
{"lineNum":"  171","line":"          Calling this function more than once will result in an exception"},
{"lineNum":"  172","line":"          being thrown."},
{"lineNum":"  173","line":""},
{"lineNum":"  174","line":"          @param args The arguments to the constructor for T"},
{"lineNum":"  175","line":"          @throw Exception If called more than once */"},
{"lineNum":"  176","line":"      template <class ... Args>"},
{"lineNum":"  177","line":"      void operator()( Args && ... args );"},
{"lineNum":"  178","line":"      // implementation deferred due to reliance on cereal::access"},
{"lineNum":"  179","line":""},
{"lineNum":"  180","line":"      //! Get a reference to the initialized underlying object"},
{"lineNum":"  181","line":"      /*! This must be called after the object has been initialized."},
{"lineNum":"  182","line":""},
{"lineNum":"  183","line":"          @return A reference to the initialized object"},
{"lineNum":"  184","line":"          @throw Exception If called before initialization */"},
{"lineNum":"  185","line":"      T * operator->()"},
{"lineNum":"  186","line":"      {"},
{"lineNum":"  187","line":"        if( !itsValid )"},
{"lineNum":"  188","line":"          throw Exception(\"Object must be initialized prior to accessing members\");"},
{"lineNum":"  189","line":""},
{"lineNum":"  190","line":"        return itsPtr;"},
{"lineNum":"  191","line":"      }"},
{"lineNum":"  192","line":""},
{"lineNum":"  193","line":"      //! Returns a raw pointer to the initialized underlying object"},
{"lineNum":"  194","line":"      /*! This is mainly intended for use with passing an instance of"},
{"lineNum":"  195","line":"          a constructed object to cereal::base_class."},
{"lineNum":"  196","line":""},
{"lineNum":"  197","line":"          It is strongly recommended to avoid using this function in"},
{"lineNum":"  198","line":"          any other circumstance."},
{"lineNum":"  199","line":""},
{"lineNum":"  200","line":"          @return A raw pointer to the initialized type */"},
{"lineNum":"  201","line":"      T * ptr()"},
{"lineNum":"  202","line":"      {"},
{"lineNum":"  203","line":"        return operator->();"},
{"lineNum":"  204","line":"      }"},
{"lineNum":"  205","line":""},
{"lineNum":"  206","line":"    private:"},
{"lineNum":"  207","line":"      template <class Ar, class TT> friend struct ::cereal::memory_detail::LoadAndConstructLoadWrapper;"},
{"lineNum":"  208","line":"      template <class Ar, class TT> friend struct ::cereal::boost_variant_detail::LoadAndConstructLoadWrapper;"},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"      construct( T * p ) : itsPtr( p ), itsEnableSharedRestoreFunction( [](){} ), itsValid( false ) {}"},
{"lineNum":"  211","line":"      construct( T * p, std::function<void()> enableSharedFunc ) : // g++4.7 ice with default lambda to std func"},
{"lineNum":"  212","line":"        itsPtr( p ), itsEnableSharedRestoreFunction( enableSharedFunc ), itsValid( false ) {}"},
{"lineNum":"  213","line":"      construct( construct const & ) = delete;"},
{"lineNum":"  214","line":"      construct & operator=( construct const & ) = delete;"},
{"lineNum":"  215","line":""},
{"lineNum":"  216","line":"      T * itsPtr;"},
{"lineNum":"  217","line":"      std::function<void()> itsEnableSharedRestoreFunction;"},
{"lineNum":"  218","line":"      bool itsValid;"},
{"lineNum":"  219","line":"  };"},
{"lineNum":"  220","line":""},
{"lineNum":"  221","line":"  // ######################################################################"},
{"lineNum":"  222","line":"  //! A class that can be made a friend to give cereal access to non public functions"},
{"lineNum":"  223","line":"  /*! If you desire non-public serialization functions within a class, cereal can only"},
{"lineNum":"  224","line":"      access these if you declare cereal::access a friend."},
{"lineNum":"  225","line":""},
{"lineNum":"  226","line":"      @code{.cpp}"},
{"lineNum":"  227","line":"      class MyClass"},
{"lineNum":"  228","line":"      {"},
{"lineNum":"  229","line":"        private:"},
{"lineNum":"  230","line":"          friend class cereal::access; // gives access to the private serialize"},
{"lineNum":"  231","line":""},
{"lineNum":"  232","line":"          template <class Archive>"},
{"lineNum":"  233","line":"          void serialize( Archive & ar )"},
{"lineNum":"  234","line":"          {"},
{"lineNum":"  235","line":"            // some code"},
{"lineNum":"  236","line":"          }"},
{"lineNum":"  237","line":"      };"},
{"lineNum":"  238","line":"      @endcode"},
{"lineNum":"  239","line":"      @ingroup Access */"},
{"lineNum":"  240","line":"  class access"},
{"lineNum":"  241","line":"  {"},
{"lineNum":"  242","line":"    public:"},
{"lineNum":"  243","line":"      // ####### Standard Serialization ########################################"},
{"lineNum":"  244","line":"      template<class Archive, class T> inline"},
{"lineNum":"  245","line":"      static auto member_serialize(Archive & ar, T & t) -> decltype(t.CEREAL_SERIALIZE_FUNCTION_NAME(ar))"},
{"lineNum":"  246","line":"      { return t.CEREAL_SERIALIZE_FUNCTION_NAME(ar); }","class":"lineNoCov","hits":"0","possible_hits":"282",},
{"lineNum":"  247","line":""},
{"lineNum":"  248","line":"      template<class Archive, class T> inline"},
{"lineNum":"  249","line":"      static auto member_save(Archive & ar, T const & t) -> decltype(t.CEREAL_SAVE_FUNCTION_NAME(ar))"},
{"lineNum":"  250","line":"      { return t.CEREAL_SAVE_FUNCTION_NAME(ar); }","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  251","line":""},
{"lineNum":"  252","line":"      template<class Archive, class T> inline"},
{"lineNum":"  253","line":"      static auto member_save_non_const(Archive & ar, T & t) -> decltype(t.CEREAL_SAVE_FUNCTION_NAME(ar))"},
{"lineNum":"  254","line":"      { return t.CEREAL_SAVE_FUNCTION_NAME(ar); }"},
{"lineNum":"  255","line":""},
{"lineNum":"  256","line":"      template<class Archive, class T> inline"},
{"lineNum":"  257","line":"      static auto member_load(Archive & ar, T & t) -> decltype(t.CEREAL_LOAD_FUNCTION_NAME(ar))"},
{"lineNum":"  258","line":"      { return t.CEREAL_LOAD_FUNCTION_NAME(ar); }","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  259","line":""},
{"lineNum":"  260","line":"      template<class Archive, class T> inline"},
{"lineNum":"  261","line":"      static auto member_save_minimal(Archive const & ar, T const & t) -> decltype(t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar))"},
{"lineNum":"  262","line":"      { return t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar); }"},
{"lineNum":"  263","line":""},
{"lineNum":"  264","line":"      template<class Archive, class T> inline"},
{"lineNum":"  265","line":"      static auto member_save_minimal_non_const(Archive const & ar, T & t) -> decltype(t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar))"},
{"lineNum":"  266","line":"      { return t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar); }"},
{"lineNum":"  267","line":""},
{"lineNum":"  268","line":"      template<class Archive, class T, class U> inline"},
{"lineNum":"  269","line":"      static auto member_load_minimal(Archive const & ar, T & t, U && u) -> decltype(t.CEREAL_LOAD_MINIMAL_FUNCTION_NAME(ar, std::forward<U>(u)))"},
{"lineNum":"  270","line":"      { return t.CEREAL_LOAD_MINIMAL_FUNCTION_NAME(ar, std::forward<U>(u)); }"},
{"lineNum":"  271","line":""},
{"lineNum":"  272","line":"      // ####### Versioned Serialization #######################################"},
{"lineNum":"  273","line":"      template<class Archive, class T> inline"},
{"lineNum":"  274","line":"      static auto member_serialize(Archive & ar, T & t, const std::uint32_t version ) -> decltype(t.CEREAL_SERIALIZE_FUNCTION_NAME(ar, version))"},
{"lineNum":"  275","line":"      { return t.CEREAL_SERIALIZE_FUNCTION_NAME(ar, version); }"},
{"lineNum":"  276","line":""},
{"lineNum":"  277","line":"      template<class Archive, class T> inline"},
{"lineNum":"  278","line":"      static auto member_save(Archive & ar, T const & t, const std::uint32_t version ) -> decltype(t.CEREAL_SAVE_FUNCTION_NAME(ar, version))"},
{"lineNum":"  279","line":"      { return t.CEREAL_SAVE_FUNCTION_NAME(ar, version); }"},
{"lineNum":"  280","line":""},
{"lineNum":"  281","line":"      template<class Archive, class T> inline"},
{"lineNum":"  282","line":"      static auto member_save_non_const(Archive & ar, T & t, const std::uint32_t version ) -> decltype(t.CEREAL_SAVE_FUNCTION_NAME(ar, version))"},
{"lineNum":"  283","line":"      { return t.CEREAL_SAVE_FUNCTION_NAME(ar, version); }"},
{"lineNum":"  284","line":""},
{"lineNum":"  285","line":"      template<class Archive, class T> inline"},
{"lineNum":"  286","line":"      static auto member_load(Archive & ar, T & t, const std::uint32_t version ) -> decltype(t.CEREAL_LOAD_FUNCTION_NAME(ar, version))"},
{"lineNum":"  287","line":"      { return t.CEREAL_LOAD_FUNCTION_NAME(ar, version); }"},
{"lineNum":"  288","line":""},
{"lineNum":"  289","line":"      template<class Archive, class T> inline"},
{"lineNum":"  290","line":"      static auto member_save_minimal(Archive const & ar, T const & t, const std::uint32_t version) -> decltype(t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar, version))"},
{"lineNum":"  291","line":"      { return t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar, version); }"},
{"lineNum":"  292","line":""},
{"lineNum":"  293","line":"      template<class Archive, class T> inline"},
{"lineNum":"  294","line":"      static auto member_save_minimal_non_const(Archive const & ar, T & t, const std::uint32_t version) -> decltype(t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar, version))"},
{"lineNum":"  295","line":"      { return t.CEREAL_SAVE_MINIMAL_FUNCTION_NAME(ar, version); }"},
{"lineNum":"  296","line":""},
{"lineNum":"  297","line":"      template<class Archive, class T, class U> inline"},
{"lineNum":"  298","line":"      static auto member_load_minimal(Archive const & ar, T & t, U && u, const std::uint32_t version) -> decltype(t.CEREAL_LOAD_MINIMAL_FUNCTION_NAME(ar, std::forward<U>(u), version))"},
{"lineNum":"  299","line":"      { return t.CEREAL_LOAD_MINIMAL_FUNCTION_NAME(ar, std::forward<U>(u), version); }"},
{"lineNum":"  300","line":""},
{"lineNum":"  301","line":"      // ####### Other Functionality ##########################################"},
{"lineNum":"  302","line":"      // for detecting inheritance from enable_shared_from_this"},
{"lineNum":"  303","line":"      template <class T> inline"},
{"lineNum":"  304","line":"      static auto shared_from_this(T & t) -> decltype(t.shared_from_this());"},
{"lineNum":"  305","line":""},
{"lineNum":"  306","line":"      // for placement new"},
{"lineNum":"  307","line":"      template <class T, class ... Args> inline"},
{"lineNum":"  308","line":"      static void construct( T *& ptr, Args && ... args )"},
{"lineNum":"  309","line":"      {"},
{"lineNum":"  310","line":"        new (ptr) T( std::forward<Args>( args )... );"},
{"lineNum":"  311","line":"      }"},
{"lineNum":"  312","line":""},
{"lineNum":"  313","line":"      // for non-placement new with a default constructor"},
{"lineNum":"  314","line":"      template <class T> inline"},
{"lineNum":"  315","line":"      static T * construct()"},
{"lineNum":"  316","line":"      {","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  317","line":"        return new T();","class":"lineNoCov","hits":"0","possible_hits":"219",},
{"lineNum":"  318","line":"      }","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  319","line":""},
{"lineNum":"  320","line":"      template <class T> inline"},
{"lineNum":"  321","line":"      static std::false_type load_and_construct(...)"},
{"lineNum":"  322","line":"      { return std::false_type(); }"},
{"lineNum":"  323","line":""},
{"lineNum":"  324","line":"      template<class T, class Archive> inline"},
{"lineNum":"  325","line":"      static auto load_and_construct(Archive & ar, ::cereal::construct<T> & construct) -> decltype(T::load_and_construct(ar, construct))"},
{"lineNum":"  326","line":"      {"},
{"lineNum":"  327","line":"        T::load_and_construct( ar, construct );"},
{"lineNum":"  328","line":"      }"},
{"lineNum":"  329","line":""},
{"lineNum":"  330","line":"      template<class T, class Archive> inline"},
{"lineNum":"  331","line":"      static auto load_and_construct(Archive & ar, ::cereal::construct<T> & construct, const std::uint32_t version) -> decltype(T::load_and_construct(ar, construct, version))"},
{"lineNum":"  332","line":"      {"},
{"lineNum":"  333","line":"        T::load_and_construct( ar, construct, version );"},
{"lineNum":"  334","line":"      }"},
{"lineNum":"  335","line":"  }; // end class access"},
{"lineNum":"  336","line":""},
{"lineNum":"  337","line":"  // ######################################################################"},
{"lineNum":"  338","line":"  // Deferred Implementation, see construct for more information"},
{"lineNum":"  339","line":"  template <class T> template <class ... Args> inline"},
{"lineNum":"  340","line":"  void construct<T>::operator()( Args && ... args )"},
{"lineNum":"  341","line":"  {"},
{"lineNum":"  342","line":"    if( itsValid )"},
{"lineNum":"  343","line":"      throw Exception(\"Attempting to construct an already initialized object\");"},
{"lineNum":"  344","line":""},
{"lineNum":"  345","line":"    ::cereal::access::construct( itsPtr, std::forward<Args>( args )... );"},
{"lineNum":"  346","line":"    itsEnableSharedRestoreFunction();"},
{"lineNum":"  347","line":"    itsValid = true;"},
{"lineNum":"  348","line":"  }"},
{"lineNum":"  349","line":"} // namespace cereal"},
{"lineNum":"  350","line":""},
{"lineNum":"  351","line":"#endif // CEREAL_ACCESS_HPP_"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 6, "covered" : 0,};
var merged_data = [];
