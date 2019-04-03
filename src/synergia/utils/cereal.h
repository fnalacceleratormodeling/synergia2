#ifndef SYNERGIA_UTILS_CEREAL_H_
#define SYNERGIA_UTILS_CEREAL_H_

#pragma clang diagnostic ignored "-Wc++11-extensions"

#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

//#include <cereal/archive/xml_iarchive.hpp>
//#include <cereal/archive/xml_oarchive.hpp>

#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/memory.hpp>

#if 0
#include <cereal/types/version.hpp>
#include <cereal/types/base_object.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/assume_abstract.hpp>
#include <cereal/types/nvp.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/shared_ptr.hpp>
#include <cereal/types/export.hpp>
#endif

#endif /* CEREAL_H */
