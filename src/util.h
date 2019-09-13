#ifndef _vgi_util_h_
#define _vgi_util_h_

#include "tuple.h"
#include <map>
#include <string>
#include <fstream>
using std::string;

namespace vgi {

struct Dictionary {
    bool load(const char * filename);
    int get_int(const char * name) const;
    double get_real(const char * name) const;
    string get_string(const char * name) const;
    template<int TDim, typename TType, Tuple_Type TClass>
    Tuple<TDim,TType,TClass> get_tuple(const char * name) const;
private:
    std::map<string, string> data;
};

} // namespace vgi

#endif // _vgi_util_h_