#ifndef _vgi_geom_h_
#define _vgi_geom_h_

#include "tuple.h"
#include <vector>

namespace vgi {
typedef unsigned int index_t;
/*
struct Material { 
    Color3d emission;
    Color3d diffuse;
    Color3d specular_color;
    double reflexion_factor;
    double refraction_factor;
    double shininess;
};
*/
struct Face {
    index_t vertices[3];
    index_t normals[3];
    index_t uvcoords[3];
    // index_t material;
};

struct Mesh {
    std::vector<Point3d> vertices;
    std::vector<Direction3d> normals;
    std::vector<Point2d> uvmap;
    // std::vector<Material> materials;
    std::vector<Face> faces;
    Triangle3d get_triangle(index_t face_index) const;
};

struct Hit {
    Point3d point;
    Direction3d normal;
    // Material material;
    Point2d uv;
    Mesh * mesh;
    index_t face_index;
};

Hit intersect(const Ray3d& ray, const Triangle3d& triangle);
Hit intersect(const Ray3d& ray, const Mesh& mesh);

} // namespace vgi
#endif // _vgi_geom_h_