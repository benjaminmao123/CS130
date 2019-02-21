#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    TODO; // determine the color
	Ray reflect;

	reflect.endpoint = intersection_point;
	reflect.direction = ray.direction - (2 * dot(ray.direction, normal) * normal);

	color = reflectivity * world.Cast_Ray(reflect, recursion_depth + 1) + (1 - reflectivity) *
			shader->Shade_Surface(ray, intersection_point, normal, recursion_depth);

    return color;
}
