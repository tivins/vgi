#include "geom.h"

namespace vgi {

bool Mesh::load(const char * filename) {
    // vertices
    // normals
    // uvmap
    // faces
    return true;
}

Triangle3d Mesh::get_triangle(index_t face_index) const {
    Triangle3d out;
    if (face_index > faces.size() - 1) return out;
    out.vertices[0] = vertices[faces[face_index][0]];
    out.vertices[1] = vertices[faces[face_index][1]];
    out.vertices[2] = vertices[faces[face_index][2]];
    return out;
}

Hit intersect(const Ray3d& ray, const Triangle3d& triangle) {
    Hit out;
    /*
    Point3d point;
    Direction3d normal;
    Point2d uv;
    Mesh * mesh;
    index_t face_index;
    */
    return out;
}

Hit intersect(const Ray3d& ray, const Mesh& mesh) {
    Hit out, tmp;
    for (index_t idx = 0; idx < mesh.faces.size(); idx++) {
        tmp = intersect(ray, mesh.get_triangle(idx));
    }
    return out;
}

} // namespace