#ifndef _vgi_scene_h_
#define _vgi_scene_h_

#include "geom.h"
#include <map>
#include <string>

namespace vgi {
using std::string;

/**
 *
 */
struct Camera {
    Point3d position, look_at;
    Direction3d up_dir;
private:
    Matrix4d proj, view;
};

/**
 *
 */
struct Material {
    Color3d emission;
    Color3d diffuse;
    Color3d specular_color;
    double reflexion_factor;
    double refraction_factor;
    double shininess;
    Image3d * color_map;
};

/**
 *
 */
struct Object {
    Object * parent = nullptr;
    Mesh * mesh = nullptr;
    Material * material = nullptr;
    inline void apply_transform(const Matrix4d& transf) {
        transform = transform * transf;
        // itransform = transform.get_inverse(); // todo
    }
private:
    Matrix4d transform, itransform;
};

/**
 *
 */
struct Scene {
    std::map<string, Object*> objects;
    std::map<string, Material*> materials;
    Camera * camera = nullptr;
    std::vector<Object*> get_emissive_objects();
    bool load(const char * filename);
};

} // namespace vgi
#endif // _vgi_scene_h_