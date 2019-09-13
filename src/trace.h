#ifndef _vgi_trace_h_
#define _vgi_trace_h_

#include "scene.h"

namespace vgi {

struct Tracer {
    Image3d output;
    Scene scene;

    inline Size2i get_paper_size() const { return output.get_size(); }
};

} // namespace vgi
#endif // _vgi_trace_h_