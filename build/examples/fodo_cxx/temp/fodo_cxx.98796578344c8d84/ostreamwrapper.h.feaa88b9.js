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
{"lineNum":"   15","line":"#ifndef CEREAL_RAPIDJSON_OSTREAMWRAPPER_H_"},
{"lineNum":"   16","line":"#define CEREAL_RAPIDJSON_OSTREAMWRAPPER_H_"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"#include \"stream.h\""},
{"lineNum":"   19","line":"#include <iosfwd>"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"#ifdef __clang__"},
{"lineNum":"   22","line":"CEREAL_RAPIDJSON_DIAG_PUSH"},
{"lineNum":"   23","line":"CEREAL_RAPIDJSON_DIAG_OFF(padded)"},
{"lineNum":"   24","line":"#endif"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"CEREAL_RAPIDJSON_NAMESPACE_BEGIN"},
{"lineNum":"   27","line":""},
{"lineNum":"   28","line":"//! Wrapper of \\c std::basic_ostream into RapidJSON\'s Stream concept."},
{"lineNum":"   29","line":"/*!"},
{"lineNum":"   30","line":"    The classes can be wrapped including but not limited to:"},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":"    - \\c std::ostringstream"},
{"lineNum":"   33","line":"    - \\c std::stringstream"},
{"lineNum":"   34","line":"    - \\c std::wpstringstream"},
{"lineNum":"   35","line":"    - \\c std::wstringstream"},
{"lineNum":"   36","line":"    - \\c std::ifstream"},
{"lineNum":"   37","line":"    - \\c std::fstream"},
{"lineNum":"   38","line":"    - \\c std::wofstream"},
{"lineNum":"   39","line":"    - \\c std::wfstream"},
{"lineNum":"   40","line":""},
{"lineNum":"   41","line":"    \\tparam StreamType Class derived from \\c std::basic_ostream."},
{"lineNum":"   42","line":"*/"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"template <typename StreamType>"},
{"lineNum":"   45","line":"class BasicOStreamWrapper {"},
{"lineNum":"   46","line":"public:"},
{"lineNum":"   47","line":"    typedef typename StreamType::char_type Ch;"},
{"lineNum":"   48","line":"    BasicOStreamWrapper(StreamType& stream) : stream_(stream) {}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"    void Put(Ch c) {"},
{"lineNum":"   51","line":"        stream_.put(c);","class":"lineNoCov","hits":"0","possible_hits":"279",},
{"lineNum":"   52","line":"    }"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":"    void Flush() {"},
{"lineNum":"   55","line":"        stream_.flush();","class":"lineNoCov","hits":"0","possible_hits":"215",},
{"lineNum":"   56","line":"    }"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"    // Not implemented"},
{"lineNum":"   59","line":"    char Peek() const { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"   60","line":"    char Take() { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"   61","line":"    size_t Tell() const { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"   62","line":"    char* PutBegin() { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"   63","line":"    size_t PutEnd(char*) { CEREAL_RAPIDJSON_ASSERT(false); return 0; }"},
{"lineNum":"   64","line":""},
{"lineNum":"   65","line":"private:"},
{"lineNum":"   66","line":"    BasicOStreamWrapper(const BasicOStreamWrapper&);"},
{"lineNum":"   67","line":"    BasicOStreamWrapper& operator=(const BasicOStreamWrapper&);"},
{"lineNum":"   68","line":""},
{"lineNum":"   69","line":"    StreamType& stream_;"},
{"lineNum":"   70","line":"};"},
{"lineNum":"   71","line":""},
{"lineNum":"   72","line":"typedef BasicOStreamWrapper<std::ostream> OStreamWrapper;"},
{"lineNum":"   73","line":"typedef BasicOStreamWrapper<std::wostream> WOStreamWrapper;"},
{"lineNum":"   74","line":""},
{"lineNum":"   75","line":"#ifdef __clang__"},
{"lineNum":"   76","line":"CEREAL_RAPIDJSON_DIAG_POP"},
{"lineNum":"   77","line":"#endif"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":"CEREAL_RAPIDJSON_NAMESPACE_END"},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"#endif // CEREAL_RAPIDJSON_OSTREAMWRAPPER_H_"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 3, "covered" : 0,};
var merged_data = [];
