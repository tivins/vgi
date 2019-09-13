#ifndef _vgi_util_h_
#define _vgi_util_h_

#include <map>
#include <string>
using std::string;

namespace vgi {

struct Dict {
    bool load(const char * filename);
    int get_int(const char * name) const;
    double get_real(const char * name) const;
    string get_string(const char * name) const;
};

} // namespace vgi

#endif // _vgi_util_h_