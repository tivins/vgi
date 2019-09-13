#include "util.h"

namespace vgi {

bool Dictionary::load(const char * filename) {
    std::ifstream file(filename);
    file.close();
}
/*
int get_int(const char * name) const;
double get_real(const char * name) const;
string get_string(const char * name) const;
template<int TDim, typename TType, Tuple_Type TClass>
Tuple<TDim,TType,TClass> get_tuple(const char * name) const;
*/

} // namespace vgi