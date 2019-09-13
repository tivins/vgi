#include "geom.h"

namespace vgi {

Triangle3d Mesh::get_triangle(index_t face_index) const {
    Triangle3d out;
    return out;
}

Hit intersect(const Ray3d& ray, const Triangle3d& triangle) {
    Hit out;
    return out;
}

Hit intersect(const Ray3d& ray, const Mesh& mesh) {
    Hit out;
    return out;
}

} // namespace