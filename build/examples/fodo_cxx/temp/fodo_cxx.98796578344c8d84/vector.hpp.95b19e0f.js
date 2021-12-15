var data = {lines:[
{"lineNum":"    1","line":"/*! \\file vector.hpp"},
{"lineNum":"    2","line":"    \\brief Support for types found in \\<vector\\>"},
{"lineNum":"    3","line":"    \\ingroup STLSupport */"},
{"lineNum":"    4","line":"/*"},
{"lineNum":"    5","line":"  Copyright (c) 2014, Randolph Voorhies, Shane Grant"},
{"lineNum":"    6","line":"  All rights reserved."},
{"lineNum":"    7","line":""},
{"lineNum":"    8","line":"  Redistribution and use in source and binary forms, with or without"},
{"lineNum":"    9","line":"  modification, are permitted provided that the following conditions are met:"},
{"lineNum":"   10","line":"      * Redistributions of source code must retain the above copyright"},
{"lineNum":"   11","line":"        notice, this list of conditions and the following disclaimer."},
{"lineNum":"   12","line":"      * Redistributions in binary form must reproduce the above copyright"},
{"lineNum":"   13","line":"        notice, this list of conditions and the following disclaimer in the"},
{"lineNum":"   14","line":"        documentation and/or other materials provided with the distribution."},
{"lineNum":"   15","line":"      * Neither the name of cereal nor the"},
{"lineNum":"   16","line":"        names of its contributors may be used to endorse or promote products"},
{"lineNum":"   17","line":"        derived from this software without specific prior written permission."},
{"lineNum":"   18","line":""},
{"lineNum":"   19","line":"  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND"},
{"lineNum":"   20","line":"  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED"},
{"lineNum":"   21","line":"  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE"},
{"lineNum":"   22","line":"  DISCLAIMED. IN NO EVENT SHALL RANDOLPH VOORHIES OR SHANE GRANT BE LIABLE FOR ANY"},
{"lineNum":"   23","line":"  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES"},
{"lineNum":"   24","line":"  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;"},
{"lineNum":"   25","line":"  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND"},
{"lineNum":"   26","line":"  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT"},
{"lineNum":"   27","line":"  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"},
{"lineNum":"   28","line":"  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."},
{"lineNum":"   29","line":"*/"},
{"lineNum":"   30","line":"#ifndef CEREAL_TYPES_VECTOR_HPP_"},
{"lineNum":"   31","line":"#define CEREAL_TYPES_VECTOR_HPP_"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"#include \"cereal/cereal.hpp\""},
{"lineNum":"   34","line":"#include <vector>"},
{"lineNum":"   35","line":""},
{"lineNum":"   36","line":"namespace cereal"},
{"lineNum":"   37","line":"{"},
{"lineNum":"   38","line":"  //! Serialization for std::vectors of arithmetic (but not bool) using binary serialization, if supported"},
{"lineNum":"   39","line":"  template <class Archive, class T, class A> inline"},
{"lineNum":"   40","line":"  typename std::enable_if<traits::is_output_serializable<BinaryData<T>, Archive>::value"},
{"lineNum":"   41","line":"                          && std::is_arithmetic<T>::value && !std::is_same<T, bool>::value, void>::type"},
{"lineNum":"   42","line":"  CEREAL_SAVE_FUNCTION_NAME( Archive & ar, std::vector<T, A> const & vector )"},
{"lineNum":"   43","line":"  {"},
{"lineNum":"   44","line":"    ar( make_size_tag( static_cast<size_type>(vector.size()) ) ); // number of elements"},
{"lineNum":"   45","line":"    ar( binary_data( vector.data(), vector.size() * sizeof(T) ) );"},
{"lineNum":"   46","line":"  }"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"  //! Serialization for std::vectors of arithmetic (but not bool) using binary serialization, if supported"},
{"lineNum":"   49","line":"  template <class Archive, class T, class A> inline"},
{"lineNum":"   50","line":"  typename std::enable_if<traits::is_input_serializable<BinaryData<T>, Archive>::value"},
{"lineNum":"   51","line":"                          && std::is_arithmetic<T>::value && !std::is_same<T, bool>::value, void>::type"},
{"lineNum":"   52","line":"  CEREAL_LOAD_FUNCTION_NAME( Archive & ar, std::vector<T, A> & vector )"},
{"lineNum":"   53","line":"  {"},
{"lineNum":"   54","line":"    size_type vectorSize;"},
{"lineNum":"   55","line":"    ar( make_size_tag( vectorSize ) );"},
{"lineNum":"   56","line":""},
{"lineNum":"   57","line":"    vector.resize( static_cast<std::size_t>( vectorSize ) );"},
{"lineNum":"   58","line":"    ar( binary_data( vector.data(), static_cast<std::size_t>( vectorSize ) * sizeof(T) ) );"},
{"lineNum":"   59","line":"  }"},
{"lineNum":"   60","line":""},
{"lineNum":"   61","line":"  //! Serialization for non-arithmetic vector types"},
{"lineNum":"   62","line":"  template <class Archive, class T, class A> inline"},
{"lineNum":"   63","line":"  typename std::enable_if<(!traits::is_output_serializable<BinaryData<T>, Archive>::value"},
{"lineNum":"   64","line":"                          || !std::is_arithmetic<T>::value) && !std::is_same<T, bool>::value, void>::type"},
{"lineNum":"   65","line":"  CEREAL_SAVE_FUNCTION_NAME( Archive & ar, std::vector<T, A> const & vector )"},
{"lineNum":"   66","line":"  {"},
{"lineNum":"   67","line":"    ar( make_size_tag( static_cast<size_type>(vector.size()) ) ); // number of elements","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   68","line":"    for(auto && v : vector)","class":"lineNoCov","hits":"0","possible_hits":"23",},
{"lineNum":"   69","line":"      ar( v );"},
{"lineNum":"   70","line":"  }"},
{"lineNum":"   71","line":""},
{"lineNum":"   72","line":"  //! Serialization for non-arithmetic vector types"},
{"lineNum":"   73","line":"  template <class Archive, class T, class A> inline"},
{"lineNum":"   74","line":"  typename std::enable_if<(!traits::is_input_serializable<BinaryData<T>, Archive>::value"},
{"lineNum":"   75","line":"                          || !std::is_arithmetic<T>::value) && !std::is_same<T, bool>::value, void>::type"},
{"lineNum":"   76","line":"  CEREAL_LOAD_FUNCTION_NAME( Archive & ar, std::vector<T, A> & vector )"},
{"lineNum":"   77","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"   78","line":"    size_type size;"},
{"lineNum":"   79","line":"    ar( make_size_tag( size ) );"},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"    vector.resize( static_cast<std::size_t>( size ) );","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"   82","line":"    for(auto && v : vector)","class":"lineNoCov","hits":"0","possible_hits":"16",},
{"lineNum":"   83","line":"      ar( v );"},
{"lineNum":"   84","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"10",},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"  //! Serialization for bool vector types"},
{"lineNum":"   87","line":"  template <class Archive, class A> inline"},
{"lineNum":"   88","line":"  void CEREAL_SAVE_FUNCTION_NAME( Archive & ar, std::vector<bool, A> const & vector )"},
{"lineNum":"   89","line":"  {"},
{"lineNum":"   90","line":"    ar( make_size_tag( static_cast<size_type>(vector.size()) ) ); // number of elements"},
{"lineNum":"   91","line":"    for(const auto v : vector)"},
{"lineNum":"   92","line":"      ar( static_cast<bool>(v) );"},
{"lineNum":"   93","line":"  }"},
{"lineNum":"   94","line":""},
{"lineNum":"   95","line":"  //! Serialization for bool vector types"},
{"lineNum":"   96","line":"  template <class Archive, class A> inline"},
{"lineNum":"   97","line":"  void CEREAL_LOAD_FUNCTION_NAME( Archive & ar, std::vector<bool, A> & vector )"},
{"lineNum":"   98","line":"  {"},
{"lineNum":"   99","line":"    size_type size;"},
{"lineNum":"  100","line":"    ar( make_size_tag( size ) );"},
{"lineNum":"  101","line":""},
{"lineNum":"  102","line":"    vector.resize( static_cast<std::size_t>( size ) );"},
{"lineNum":"  103","line":"    for(auto v : vector)"},
{"lineNum":"  104","line":"    {"},
{"lineNum":"  105","line":"      bool b;"},
{"lineNum":"  106","line":"      ar( b );"},
{"lineNum":"  107","line":"      v = b;"},
{"lineNum":"  108","line":"    }"},
{"lineNum":"  109","line":"  }"},
{"lineNum":"  110","line":"} // namespace cereal"},
{"lineNum":"  111","line":""},
{"lineNum":"  112","line":"#endif // CEREAL_TYPES_VECTOR_HPP_"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 6, "covered" : 0,};
var merged_data = [];
