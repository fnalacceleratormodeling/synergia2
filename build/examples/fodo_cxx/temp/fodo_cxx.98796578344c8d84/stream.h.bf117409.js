var data = {lines:[
{"lineNum":"    1","line":"// Tencent is pleased to support the open source community by making RapidJSON available."},
{"lineNum":"    2","line":"//"},
{"lineNum":"    3","line":"// Copyright (C) 2015 THL A29 Limited, a Tencent company, and Milo Yip. All rights reserved."},
{"lineNum":"    4","line":"//"},
{"lineNum":"    5","line":"// Licensed under the MIT License (the \"License\"); you may not use this file except"},
{"lineNum":"    6","line":"// in compliance with the License. You may obtain a copy of the License at"},
{"lineNum":"    7","line":"//"},
{"lineNum":"    8","line":"// http://opensource.org/licenses/MIT"},
{"lineNum":"    9","line":"//"},
{"lineNum":"   10","line":"// Unless required by applicable law or agreed to in writing, software distributed"},
{"lineNum":"   11","line":"// under the License is distributed on an \"AS IS\" BASIS, WITHOUT WARRANTIES OR"},
{"lineNum":"   12","line":"// CONDITIONS OF ANY KIND, either express or implied. See the License for the"},
{"lineNum":"   13","line":"// specific language governing permissions and limitations under the License."},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"#include \"rapidjson.h\""},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"#ifndef CEREAL_RAPIDJSON_STREAM_H_"},
{"lineNum":"   18","line":"#define CEREAL_RAPIDJSON_STREAM_H_"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"#include \"encodings.h\""},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"CEREAL_RAPIDJSON_NAMESPACE_BEGIN"},
{"lineNum":"   23","line":""},
{"lineNum":"   24","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   25","line":"//  Stream"},
{"lineNum":"   26","line":""},
{"lineNum":"   27","line":"/*! \\class rapidjson::Stream"},
{"lineNum":"   28","line":"    \\brief Concept for reading and writing characters."},
{"lineNum":"   29","line":""},
{"lineNum":"   30","line":"    For read-only stream, no need to implement PutBegin(), Put(), Flush() and PutEnd()."},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":"    For write-only stream, only need to implement Put() and Flush()."},
{"lineNum":"   33","line":""},
{"lineNum":"   34","line":"\\code"},
{"lineNum":"   35","line":"concept Stream {"},
{"lineNum":"   36","line":"    typename Ch;    //!< Character type of the stream."},
{"lineNum":"   37","line":""},
{"lineNum":"   38","line":"    //! Read the current character from stream without moving the read cursor."},
{"lineNum":"   39","line":"    Ch Peek() const;"},
{"lineNum":"   40","line":""},
{"lineNum":"   41","line":"    //! Read the current character from stream and moving the read cursor to next character."},
{"lineNum":"   42","line":"    Ch Take();"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"    //! Get the current read cursor."},
{"lineNum":"   45","line":"    //! \\return Number of characters read from start."},
{"lineNum":"   46","line":"    size_t Tell();"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"    //! Begin writing operation at the current read pointer."},
{"lineNum":"   49","line":"    //! \\return The begin writer pointer."},
{"lineNum":"   50","line":"    Ch* PutBegin();"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"    //! Write a character."},
{"lineNum":"   53","line":"    void Put(Ch c);"},
{"lineNum":"   54","line":""},
{"lineNum":"   55","line":"    //! Flush the buffer."},
{"lineNum":"   56","line":"    void Flush();"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"    //! End the writing operation."},
{"lineNum":"   59","line":"    //! \\param begin The begin write pointer returned by PutBegin()."},
{"lineNum":"   60","line":"    //! \\return Number of characters written."},
{"lineNum":"   61","line":"    size_t PutEnd(Ch* begin);"},
{"lineNum":"   62","line":"}"},
{"lineNum":"   63","line":"\\endcode"},
{"lineNum":"   64","line":"*/"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"//! Provides additional information for stream."},
{"lineNum":"   67","line":"/*!"},
{"lineNum":"   68","line":"    By using traits pattern, this type provides a default configuration for stream."},
{"lineNum":"   69","line":"    For custom stream, this type can be specialized for other configuration."},
{"lineNum":"   70","line":"    See TEST(Reader, CustomStringStream) in readertest.cpp for example."},
{"lineNum":"   71","line":"*/"},
{"lineNum":"   72","line":"template<typename Stream>"},
{"lineNum":"   73","line":"struct StreamTraits {"},
{"lineNum":"   74","line":"    //! Whether to make local copy of stream for optimization during parsing."},
{"lineNum":"   75","line":"    /*!"},
{"lineNum":"   76","line":"        By default, for safety, streams do not use local copy optimization."},
{"lineNum":"   77","line":"        Stream that can be copied fast should specialize this, like StreamTraits<StringStream>."},
{"lineNum":"   78","line":"    */"},
{"lineNum":"   79","line":"    enum { copyOptimization = 0 };"},
{"lineNum":"   80","line":"};"},
{"lineNum":"   81","line":""},
{"lineNum":"   82","line":"//! Reserve n characters for writing to a stream."},
{"lineNum":"   83","line":"template<typename Stream>"},
{"lineNum":"   84","line":"inline void PutReserve(Stream& stream, size_t count) {"},
{"lineNum":"   85","line":"    (void)stream;"},
{"lineNum":"   86","line":"    (void)count;"},
{"lineNum":"   87","line":"}"},
{"lineNum":"   88","line":""},
{"lineNum":"   89","line":"//! Write character to a stream, presuming buffer is reserved."},
{"lineNum":"   90","line":"template<typename Stream>"},
{"lineNum":"   91","line":"inline void PutUnsafe(Stream& stream, typename Stream::Ch c) {"},
{"lineNum":"   92","line":"    stream.Put(c);"},
{"lineNum":"   93","line":"}"},
{"lineNum":"   94","line":""},
{"lineNum":"   95","line":"//! Put N copies of a character to a stream."},
{"lineNum":"   96","line":"template<typename Stream, typename Ch>"},
{"lineNum":"   97","line":"inline void PutN(Stream& stream, Ch c, size_t n) {"},
{"lineNum":"   98","line":"    PutReserve(stream, n);"},
{"lineNum":"   99","line":"    for (size_t i = 0; i < n; i++)","class":"lineNoCov","hits":"0","possible_hits":"32",},
{"lineNum":"  100","line":"        PutUnsafe(stream, c);"},
{"lineNum":"  101","line":"}"},
{"lineNum":"  102","line":""},
{"lineNum":"  103","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  104","line":"// GenericStreamWrapper"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"//! A Stream Wrapper"},
{"lineNum":"  107","line":"/*! \\tThis string stream is a wrapper for any stream by just forwarding any"},
{"lineNum":"  108","line":"    \\treceived message to the origin stream."},
{"lineNum":"  109","line":"    \\note implements Stream concept"},
{"lineNum":"  110","line":"*/"},
{"lineNum":"  111","line":""},
{"lineNum":"  112","line":"#if defined(_MSC_VER) && _MSC_VER <= 1800"},
{"lineNum":"  113","line":"CEREAL_RAPIDJSON_DIAG_PUSH"},
{"lineNum":"  114","line":"CEREAL_RAPIDJSON_DIAG_OFF(4702)  // unreachable code"},
{"lineNum":"  115","line":"CEREAL_RAPIDJSON_DIAG_OFF(4512)  // assignment operator could not be generated"},
{"lineNum":"  116","line":"#endif"},
{"lineNum":"  117","line":""},
{"lineNum":"  118","line":"template <typename InputStream, typename Encoding = UTF8<> >"},
{"lineNum":"  119","line":"class GenericStreamWrapper {"},
{"lineNum":"  120","line":"public:"},
{"lineNum":"  121","line":"    typedef typename Encoding::Ch Ch;"},
{"lineNum":"  122","line":"    GenericStreamWrapper(InputStream& is): is_(is) {}"},
{"lineNum":"  123","line":""},
{"lineNum":"  124","line":"    Ch Peek() const { return is_.Peek(); }"},
{"lineNum":"  125","line":"    Ch Take() { return is_.Take(); }"},
{"lineNum":"  126","line":"    size_t Tell() { return is_.Tell(); }"},
{"lineNum":"  127","line":"    Ch* PutBegin() { return is_.PutBegin(); }"},
{"lineNum":"  128","line":"    void Put(Ch ch) { is_.Put(ch); }"},
{"lineNum":"  129","line":"    void Flush() { is_.Flush(); }"},
{"lineNum":"  130","line":"    size_t PutEnd(Ch* ch) { return is_.PutEnd(ch); }"},
{"lineNum":"  131","line":""},
{"lineNum":"  132","line":"    // wrapper for MemoryStream"},
{"lineNum":"  133","line":"    const Ch* Peek4() const { return is_.Peek4(); }"},
{"lineNum":"  134","line":""},
{"lineNum":"  135","line":"    // wrapper for AutoUTFInputStream"},
{"lineNum":"  136","line":"    UTFType GetType() const { return is_.GetType(); }"},
{"lineNum":"  137","line":"    bool HasBOM() const { return is_.HasBOM(); }"},
{"lineNum":"  138","line":""},
{"lineNum":"  139","line":"protected:"},
{"lineNum":"  140","line":"    InputStream& is_;"},
{"lineNum":"  141","line":"};"},
{"lineNum":"  142","line":""},
{"lineNum":"  143","line":"#if defined(_MSC_VER) && _MSC_VER <= 1800"},
{"lineNum":"  144","line":"CEREAL_RAPIDJSON_DIAG_POP"},
{"lineNum":"  145","line":"#endif"},
{"lineNum":"  146","line":""},
{"lineNum":"  147","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  148","line":"// StringStream"},
{"lineNum":"  149","line":""},
{"lineNum":"  150","line":"//! Read-only string stream."},
{"lineNum":"  151","line":"/*! \\note implements Stream concept"},
{"lineNum":"  152","line":"*/"},
{"lineNum":"  153","line":"template <typename Encoding>"},
{"lineNum":"  154","line":"struct GenericStringStream {"},
{"lineNum":"  155","line":"    typedef typename Encoding::Ch Ch;"},
{"lineNum":"  156","line":""},
{"lineNum":"  157","line":"    GenericStringStream(const Ch *src) : src_(src), head_(src) {}"},
{"lineNum":"  158","line":""},
{"lineNum":"  159","line":"    Ch Peek() const { return *src_; }","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  160","line":"    Ch Take() { return *src_++; }"},
{"lineNum":"  161","line":"    size_t Tell() const { return static_cast<size_t>(src_ - head_); }"},
{"lineNum":"  162","line":""},
{"lineNum":"  163","line":"    Ch* PutBegin() { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"  164","line":"    void Put(Ch) { CEREAL_RAPIDJSON_ASSERT(false); }"},
{"lineNum":"  165","line":"    void Flush() { CEREAL_RAPIDJSON_ASSERT(false); }"},
{"lineNum":"  166","line":"    size_t PutEnd(Ch*) { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"  167","line":""},
{"lineNum":"  168","line":"    const Ch* src_;     //!< Current read position."},
{"lineNum":"  169","line":"    const Ch* head_;    //!< Original head of the string."},
{"lineNum":"  170","line":"};"},
{"lineNum":"  171","line":""},
{"lineNum":"  172","line":"template <typename Encoding>"},
{"lineNum":"  173","line":"struct StreamTraits<GenericStringStream<Encoding> > {"},
{"lineNum":"  174","line":"    enum { copyOptimization = 1 };"},
{"lineNum":"  175","line":"};"},
{"lineNum":"  176","line":""},
{"lineNum":"  177","line":"//! String stream with UTF8 encoding."},
{"lineNum":"  178","line":"typedef GenericStringStream<UTF8<> > StringStream;"},
{"lineNum":"  179","line":""},
{"lineNum":"  180","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  181","line":"// InsituStringStream"},
{"lineNum":"  182","line":""},
{"lineNum":"  183","line":"//! A read-write string stream."},
{"lineNum":"  184","line":"/*! This string stream is particularly designed for in-situ parsing."},
{"lineNum":"  185","line":"    \\note implements Stream concept"},
{"lineNum":"  186","line":"*/"},
{"lineNum":"  187","line":"template <typename Encoding>"},
{"lineNum":"  188","line":"struct GenericInsituStringStream {"},
{"lineNum":"  189","line":"    typedef typename Encoding::Ch Ch;"},
{"lineNum":"  190","line":""},
{"lineNum":"  191","line":"    GenericInsituStringStream(Ch *src) : src_(src), dst_(0), head_(src) {}"},
{"lineNum":"  192","line":""},
{"lineNum":"  193","line":"    // Read"},
{"lineNum":"  194","line":"    Ch Peek() { return *src_; }"},
{"lineNum":"  195","line":"    Ch Take() { return *src_++; }"},
{"lineNum":"  196","line":"    size_t Tell() { return static_cast<size_t>(src_ - head_); }"},
{"lineNum":"  197","line":""},
{"lineNum":"  198","line":"    // Write"},
{"lineNum":"  199","line":"    void Put(Ch c) { CEREAL_RAPIDJSON_ASSERT(dst_ != 0); *dst_++ = c; }"},
{"lineNum":"  200","line":""},
{"lineNum":"  201","line":"    Ch* PutBegin() { return dst_ = src_; }"},
{"lineNum":"  202","line":"    size_t PutEnd(Ch* begin) { return static_cast<size_t>(dst_ - begin); }"},
{"lineNum":"  203","line":"    void Flush() {}"},
{"lineNum":"  204","line":""},
{"lineNum":"  205","line":"    Ch* Push(size_t count) { Ch* begin = dst_; dst_ += count; return begin; }"},
{"lineNum":"  206","line":"    void Pop(size_t count) { dst_ -= count; }"},
{"lineNum":"  207","line":""},
{"lineNum":"  208","line":"    Ch* src_;"},
{"lineNum":"  209","line":"    Ch* dst_;"},
{"lineNum":"  210","line":"    Ch* head_;"},
{"lineNum":"  211","line":"};"},
{"lineNum":"  212","line":""},
{"lineNum":"  213","line":"template <typename Encoding>"},
{"lineNum":"  214","line":"struct StreamTraits<GenericInsituStringStream<Encoding> > {"},
{"lineNum":"  215","line":"    enum { copyOptimization = 1 };"},
{"lineNum":"  216","line":"};"},
{"lineNum":"  217","line":""},
{"lineNum":"  218","line":"//! Insitu string stream with UTF8 encoding."},
{"lineNum":"  219","line":"typedef GenericInsituStringStream<UTF8<> > InsituStringStream;"},
{"lineNum":"  220","line":""},
{"lineNum":"  221","line":"CEREAL_RAPIDJSON_NAMESPACE_END"},
{"lineNum":"  222","line":""},
{"lineNum":"  223","line":"#endif // CEREAL_RAPIDJSON_STREAM_H_"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 2, "covered" : 0,};
var merged_data = [];
